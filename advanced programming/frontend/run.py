import io
import os
import base64

import matplotlib
import pandas as pd
from werkzeug.utils import secure_filename
from flask import Flask, flash, request, redirect, url_for, render_template, session

from plab2.utils import ApiInterface
from plab2.startup import DATA_DIR
from plab2.constants import UNIPROT, HGNC
from plab2.network import Analyzer, Network, Statistics

matplotlib.use('Agg')

UPLOAD_FOLDER = os.path.join(DATA_DIR, 'uploads')
ALLOWED_EXTENSIONS = {"tsv", "csv"}

app = Flask(__name__)
app.secret_key = "someSecretKey"

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_PATH'] = 16 * 1024 * 1024


@app.route("/", methods=['GET', 'POST'])
def home():
    if request.method == 'POST':  # User choosing file type
        file_type = request.form['graph_input']
        import_check = toggle_file_type(file_type)

        if session.get('files') is not None:  # File(s) imported
            sum_stats_html = create_stats()
            return render_template('template.html', sum_stats=sum_stats_html, import_msg=import_check)

        else:
            return render_template('template.html', import_msg=import_check)

    # if session.get('files') is not None:  # File(s) already imported
    #     sum_stats_html = create_stats()
    #     return render_template('template.html', sum_stats=sum_stats_html)

    return render_template('template.html')


@app.route('/plot.png', methods=['GET', 'POST'])
def plot_png():
    analyzer = create_graph_object(object_type='analyzer')

    if isinstance(analyzer, str):  # No files found error message
        return render_template('template.html', header_msg=analyzer)

    tmp_path = os.path.join(DATA_DIR, "graph.png")
    plt = analyzer.generate_graph_image(graph_output_path=tmp_path)
    img = io.BytesIO()
    plt.savefig(img, format='png')
    plt.close()
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode('utf8')
    return render_template('plot.html', plot_url=plot_url)


def toggle_file_type(file_type: str) -> str:
    ext = ".csv" if file_type == "ppi" else ".tsv"
    files = [
        os.path.join(app.config['UPLOAD_FOLDER'], f)
        for f in os.listdir(app.config['UPLOAD_FOLDER'])
        if f.endswith(ext)
    ]

    if file_type == "ppi" and len(files) != 1:
        return "Error importing PPI file. Please clear uploaded contents and upload a single PPI file."

    elif file_type == "ne_list" and len(files) != 2:
        return "Error importing node/edge lists. Please clear uploaded contents and upload a single node list \
        and a single edge list."

    else:
        session['file_type'] = file_type
        session['files'] = files
        return "File(s) imported successfully"


@app.route("/upload")
def upload():
    return render_template('upload.html')


def create_graph_object(object_type: str):
    if 'files' not in session or len(session['files']) == 0:
        return "No graph files Found!"

    # Check what class to intialize
    if object_type == "analyzer":
        network_class = Analyzer

    elif object_type == "network":
        network_class = Network

    elif object_type == "statistics":
        network_class = Statistics

    else:
        raise ValueError("Invalid object_type passed!")

    # Check which file types to import
    if session['file_type'] == "ppi":
        ppi_file = session['files'][0]
        go = network_class(ppi_file=ppi_file)
        go.import_ppi()

    elif session['file_type'] == "ne_list":
        ne_lists = session['files']

        if ne_lists[0].lower().startswith("node"):
            node_list, edge_list = ne_lists

        else:
            edge_list, node_list = ne_lists

        go = network_class(node_list=node_list, edge_list=edge_list)

    else:
        raise ValueError(f"{session['file_type']} is not a valid file type!")

    return go


def create_stats() -> str:
    stats = create_graph_object(object_type='statistics')

    if isinstance(stats, str):  # No files found error message
        return render_template('template.html', header_msg=stats)

    stats.summary_statistics()
    stats.sum_stats.drop(5, inplace=True)  # 5 = Average Degree Connectivity
    stats_table_html = f"<h4>Summary Statistics for Imported Graph</h4>" + \
                       stats.sum_stats.to_html(header="true", table_id="table", index=False)

    return stats_table_html


@app.route("/paths", methods=['POST'])
def find_shortest_paths():
    values = [x.strip() for x in request.form['textbox'].split(",")]
    if len(values) != 2:
        error_msg_html = "More than 2 HGNC symbols were passed!"
        return render_template('template.html', paths_header=error_msg_html, sum_stats=create_stats())

    source, target = values
    analyzer = create_graph_object(object_type='analyzer')

    if isinstance(analyzer, str):  # No files found error message
        return render_template('template.html', header_msg=analyzer)

    sp = analyzer.shortest_paths(source.upper(), target.upper())

    if not sp:
        error_msg_html = f"No paths found for {source.upper()} and {target.upper()}"
        return render_template('template.html', paths_header=error_msg_html, sum_stats=create_stats())

    else:
        sp_names = [analyzer.nodes[node_id]['symbol'] for node_id in sp[0]]
        pos_html = f"Shortest path for {source.upper()} and {target.upper()}"
        return render_template('template.html', paths_header=pos_html, shortest_path=sp_names,
                               sum_stats=create_stats())


@app.route("/get_info", methods=['POST'])
def get_hgnc_info():
    """Retrieves identifier information for a given HGNC symbol."""
    symbol = request.form['textbox']
    ids = ApiInterface(symbol.upper()).get_identifiers()
    if ids:
        entries_to_display = [{
            "HGNC ID": ids[HGNC]['hgnc_id'],
            "Emsembl ID": ids[HGNC]['ensembl'],
            "UniProt ID": ids[HGNC][UNIPROT][0],
            "Tax ID": ids[UNIPROT]['tax_id'],
            "Protein Name": ids[UNIPROT]["full_name"]
        }]

        id_table = pd.DataFrame(entries_to_display)
        id_table_html = f"<h4>Results for: <b>{symbol.upper()}</b></h4>" + \
                        id_table.to_html(header="true", table_id="table", index=False)

    else:
        id_table_html = f"No information found for {symbol}"

    return render_template('template.html', symbol_name=symbol.upper(), hgnc_info=id_table_html)


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/uploader', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        if request.form['uploader_action'] == "Upload":
            # check if the post request has the file part
            if 'file' not in request.files:
                flash('No file part')
                return redirect(request.url)
            file = request.files['file']
            # if user does not select file, browser also
            # submit an empty part without filename
            if file.filename == '':
                flash('No file selected')
                return redirect(request.url)
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                return render_template('template.html', header_msg="File uploaded successfully")

            return redirect(url_for('home'))

        else:  # Clear everything
            for file in os.listdir(app.config['UPLOAD_FOLDER']):
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], file)
                os.remove(file_path)
            return render_template('upload.html', clear_msg="Uploaded files successfully removed.")


if __name__ == '__main__':
    app.run(debug=True)
