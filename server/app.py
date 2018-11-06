import os
import traceback
import uuid

import wtforms
from celery import Celery
from flask import Flask, render_template, request, Blueprint, redirect, url_for, send_file, abort
from werkzeug.utils import secure_filename
from wtforms.validators import DataRequired

from jacks.jacks_io import loadJacksFullResultsFromPickle, getSortedGenes, getGeneWs, \
    REP_HDR_DEFAULT, SAMPLE_HDR_DEFAULT, SGRNA_HDR_DEFAULT, GENE_HDR_DEFAULT, preprocess, load_data_and_run, \
    PICKLE_FILENAME, GENE_FILENAME, GRNA_FILENAME
from plot_heatmap import plot_heatmap

APP_ROOT = os.path.dirname(os.path.dirname(__file__))
ANALYSIS_FOLDER = os.path.join(APP_ROOT, 'results')
CELERY_BROKER_URL = 'CELERY_BROKER_URL'
CELERY_RESULT_BACKEND = 'CELERY_RESULT_BACKEND'
MB = 1024 * 1024

bp = Blueprint('jacks', __name__)
app = Flask(__name__, template_folder="templates", static_url_path="/JACKS/static")
app.config['MAX_CONTENT_LENGTH'] = 400 * MB
app.config[CELERY_BROKER_URL] = os.getenv(CELERY_BROKER_URL, 'redis://localhost:6379/0')
app.config[CELERY_RESULT_BACKEND] = os.getenv(CELERY_RESULT_BACKEND, 'redis://localhost:6379/0')

celery = Celery(app.name, broker=app.config[CELERY_BROKER_URL], backend=app.config[CELERY_RESULT_BACKEND])
celery.conf.update(app.config)

YUSAV1_0 = 'yusav1.0'
# YUSAV1_1 = 'yusav1.1'
AVANA = 'avana'
# BRUNELLO = 'brunello'
# WHITEHEAD = 'whitehead'
# TORONTO = 'toronto'
REFERENCE_LIB_DICT = {
    YUSAV1_0: os.path.join(APP_ROOT, 'reference_grna_efficacies', 'yusa_v10_grna_JACKS_results.txt'),
    AVANA: os.path.join(APP_ROOT, 'reference_grna_efficacies', 'avana_grna_JACKS_results.txt'),
}
ANALYSIS_FILE_DICT = {
    'gene': GENE_FILENAME,
    'grna': GRNA_FILENAME
}


@celery.task()
def run_analysis(sample_spec, gene_spec, ctrl_spec, reffile, x_ref, outprefix):
    load_data_and_run(sample_spec, gene_spec, ctrl_spec, reffile, x_ref, outprefix)


def run_jacks_async(countfile, replicatefile, guidemappingfile,
                  rep_hdr, sample_hdr, common_ctrl_sample,
                  sgrna_hdr, gene_hdr,
                  outprefix, reffile):
    sample_spec, ctrl_spec, gene_spec, x_ref = preprocess(countfile, replicatefile, guidemappingfile,
             rep_hdr, sample_hdr, common_ctrl_sample,
             sgrna_hdr=sgrna_hdr, gene_hdr=gene_hdr,
             outprefix=outprefix, reffile=reffile)
    run_analysis.delay(sample_spec, gene_spec, ctrl_spec, reffile, x_ref, outprefix)


def get_analysis_file(analysis_id, analysis_file):
    return os.path.join(ANALYSIS_FOLDER, analysis_id, analysis_file)


class JacksForm(wtforms.Form):
    raw_count_file = wtforms.FileField('Raw count file', validators=[DataRequired()])
    replicate_map_file = wtforms.FileField('Replicate map file', validators=[DataRequired()])
    header_replicates = wtforms.StringField('Header for replicates', default=REP_HDR_DEFAULT)
    header_sample = wtforms.StringField('Header for sample', default=SAMPLE_HDR_DEFAULT)
    ctrl_sample_name = wtforms.StringField('Name for a control sample', default="CTRL")
    grna_gene_map_file = wtforms.FileField('gRNA-Gene map file', validators=[DataRequired()])
    header_grna = wtforms.StringField('Header for gRNA id', default=SGRNA_HDR_DEFAULT)
    header_gene = wtforms.StringField('Header for Gene', default=GENE_HDR_DEFAULT)
    use_reference_lib = wtforms.BooleanField('Use reference library', default=False)
    reference_lib = wtforms.SelectField('Reference library',
                                        choices=[("", 'None'),
                                                 (YUSAV1_0, 'Yusa v1.0'),
                                                 # ('yusav1.1', 'Yusa v1.1'),
                                                 (AVANA, 'Avana')],
                                                 # ('brunello', 'Brunello'),
                                                 # ('whitehead', 'Whitehead'),
                                                 # ('toronto', 'Toronto Knockout')],
                                                 default="")


@app.route('/')
def hello():
    return redirect(url_for('jacks.start_analysis'))


@bp.route('/', methods=["GET", "POST"])
def start_analysis():
    template = "index.html"
    form = JacksForm(request.form)
    if request.method == 'POST':
        analysis_id = str(uuid.uuid4()).replace("-", "")[:12] + "/"

        def get_file(file):
            return os.path.join(ANALYSIS_FOLDER, analysis_id, secure_filename(request.files[file].filename))

        raw_count_file = get_file("raw_count_file")
        replicate_map_file = get_file("replicate_map_file")
        header_replicates = form.header_replicates.data
        header_sample = form.header_sample.data
        common_ctrl_sample = form.ctrl_sample_name.data
        grna_gene_map_file = get_file("grna_gene_map_file")
        header_grna = form.header_grna.data
        header_gene = form.header_gene.data
        reference_lib = form.reference_lib.data
        if reference_lib:
            reference_lib = REFERENCE_LIB_DICT[reference_lib]
        files = ['raw_count_file', 'replicate_map_file', 'grna_gene_map_file']
        for file in files:
            f = request.files[file]
            file_path = get_file(file)
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            f.save(file_path)
        try:
            run_jacks_async(countfile=raw_count_file, replicatefile=replicate_map_file,
                                guidemappingfile=grna_gene_map_file,
                                rep_hdr=header_replicates, sample_hdr=header_sample, common_ctrl_sample=common_ctrl_sample,
                                sgrna_hdr=header_grna, gene_hdr=header_gene,
                                outprefix=os.path.join(ANALYSIS_FOLDER, analysis_id), reffile=reference_lib)
        except Exception as e:
            return render_template(template, form=form, error_message=traceback.format_exc().splitlines()[-1])
        return render_template(template, form=form, analysis_id=analysis_id)
    return render_template(template, form=form)


@bp.route('/results/<path:analysis_id>', methods=["GET"])
def retrieve_results(analysis_id):
    template = "results.html"
    picklefile = get_analysis_file(analysis_id, PICKLE_FILENAME)
    if os.path.isfile(picklefile):
        jacks_results, cell_lines, gene_grnas = loadJacksFullResultsFromPickle(picklefile)
        table = []
        sorted_genes = getSortedGenes(jacks_results)
        for gene in sorted_genes:
            row = [gene[1], gene[0]]
            row.extend(getGeneWs(jacks_results, gene[1]))
            table.append(row)
        return render_template(template, table=table, cell_lines=cell_lines)
    else:
        return render_template(template)


@bp.route('/results/<analysis_id>/gene/<gene>', methods=["GET"])
def plot_gene_heatmap(analysis_id, gene):
    template = "plot.html"
    picklefile = get_analysis_file(analysis_id, PICKLE_FILENAME)
    if os.path.isfile(picklefile):
        image_path = os.path.join("results", analysis_id, "figure.png")
        full_image_path = os.path.join(os.path.dirname(__file__), "static", image_path)
        if not os.path.exists(full_image_path):
            plot_heatmap(picklefile, gene, full_image_path)
        return render_template(template, gene=gene, image_path=url_for('static', filename=image_path))
    else:
        return render_template(template, gene=gene)


@bp.route('/results/<analysis_id>/download/<file>', methods=["GET"])
def download(analysis_id, file):
    try:
        filename = get_analysis_file(analysis_id, ANALYSIS_FILE_DICT[file])
        return send_file(filename, as_attachment=True)
    except KeyError:
        return "Unknown analysis file requested", 400


app.register_blueprint(bp, url_prefix='/JACKS')

if __name__ == '__main__':
    app.run()
