import os
import uuid

import wtforms
from celery import Celery
from flask import Flask, render_template, request, Blueprint, redirect, url_for

from jacks.jacks_io import loadJacksFullResultsFromPickle, getSortedGenes, getGeneWs, runJACKS, PICKLE_FILENAME, \
    REP_HDR_DEFAULT, SAMPLE_HDR_DEFAULT, SGRNA_HDR_DEFAULT
from plot_heatmap import plot_heatmap

APP_ROOT = os.path.join(os.path.dirname(__file__), "..")
CELERY_BROKER_URL = 'CELERY_BROKER_URL'
CELERY_RESULT_BACKEND = 'CELERY_RESULT_BACKEND'

bp = Blueprint('jacks', __name__)
app = Flask(__name__, template_folder="templates", static_url_path="/JACKS/static")
app.config[CELERY_BROKER_URL] = os.getenv(CELERY_BROKER_URL, 'redis://localhost:6379/0')
app.config[CELERY_RESULT_BACKEND] = os.getenv(CELERY_RESULT_BACKEND, 'redis://localhost:6379/0')

celery = Celery(app.name, broker=app.config[CELERY_BROKER_URL], backend=app.config[CELERY_RESULT_BACKEND])
celery.conf.update(app.config)


@celery.task()
def send_to_jacks(countfile, replicatefile, guidemappingfile,
                  rep_hdr, sample_hdr, common_ctrl_sample,
                  sgrna_hdr, gene_hdr,
                  outprefix, reffile):
    runJACKS(countfile, replicatefile, guidemappingfile,
             rep_hdr, sample_hdr, common_ctrl_sample,
             sgrna_hdr=sgrna_hdr, gene_hdr=gene_hdr,
             outprefix=outprefix, reffile=reffile)


def get_pickle_file(analysis_id):
    return os.path.join(APP_ROOT, "results", analysis_id, PICKLE_FILENAME)


class JacksForm(wtforms.Form):
    raw_count_file = wtforms.FileField('Raw count file')
    replicate_map_file = wtforms.FileField('Replicate map field')
    header_replicates = wtforms.StringField('Header for replicates')
    header_sample = wtforms.StringField('Header for sample')
    ctrl_sample_name = wtforms.StringField('Name for a control sample')
    grna_gene_map_file = wtforms.FileField('gRNA-Gene map file')
    header_grna = wtforms.StringField('Header for gRNA id')
    header_gene = wtforms.StringField('Header for Gene')
    use_reference_lib = wtforms.RadioField('Use reference library', choices=[('ref', 'Use reference library'),
                                                                             ('none', 'No reference')])
    reference_lib = wtforms.SelectField('Reference library', choices=[("", 'None'),
                                                                      ('yusav1.0', 'Yusa v1.0'),
                                                                      ('yusav1.1', 'Yusa v1.1'),
                                                                      ('avana', 'Avana'),
                                                                      ('brunello', 'Brunello'),
                                                                      ('whitehead', 'Whitehead'),
                                                                      ('toronto', 'Toronto Knockout')])
    max_genes_display = wtforms.IntegerField('Max genes to display', default=20)


@app.route('/')
def hello():
    return redirect(url_for('jacks.start_analysis'))


@bp.route('/', methods=["GET", "POST"])
def start_analysis():
    template = "index.html"
    form = JacksForm(request.form)
    if request.method == 'POST':
        raw_count_file = form.raw_count_file.data
        replicate_map_file = form.replicate_map_file.data
        header_replicates = form.header_replicates.data
        header_sample = form.header_sample.data
        ctrl_sample_name = form.ctrl_sample_name.data
        grna_gene_map_file = form.grna_gene_map_file.data
        header_grna = form.header_grna.data
        header_gene = form.header_gene.data
        reference_lib = form.reference_lib.data
        if not raw_count_file:
            raw_count_file = grna_gene_map_file = os.path.join(APP_ROOT, "jacks/example-small/example_count_data.tab")
            replicate_map_file = os.path.join(APP_ROOT, "jacks/example-small/example_repmap.tab")
            reference_lib = None
            header_grna = SGRNA_HDR_DEFAULT
            header_gene = 'gene'
            header_sample = SAMPLE_HDR_DEFAULT
            header_replicates = REP_HDR_DEFAULT
            common_ctrl_sample = "CTRL"

        analysis_id = str(uuid.uuid4()).replace("-", "")[:12] + "/"
        send_to_jacks.delay(countfile=raw_count_file, replicatefile=replicate_map_file,
                            guidemappingfile=grna_gene_map_file,
                            rep_hdr=header_replicates, sample_hdr=header_sample, common_ctrl_sample=common_ctrl_sample,
                            sgrna_hdr=header_grna, gene_hdr=header_gene,
                            outprefix=os.path.join(APP_ROOT, "results/") + analysis_id, reffile=reference_lib)
        return render_template(template, form=form, analysis_id=analysis_id)
    return render_template(template, form=form)


@bp.route('/results/<path:analysis_id>', methods=["GET"])
def retrieve_results(analysis_id):
    template = "results.html"
    picklefile = get_pickle_file(analysis_id)
    if os.path.isfile(picklefile):
        jacks_results, cell_lines, gene_grnas = loadJacksFullResultsFromPickle(get_pickle_file(analysis_id))
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
    picklefile = get_pickle_file(analysis_id)
    if os.path.isfile(picklefile):
        image_path = os.path.join("results", analysis_id, "figure.png")
        full_image_path = os.path.join(APP_ROOT, "server", "static", image_path)
        plot_heatmap(picklefile, gene, full_image_path)
        return render_template(template, image_path=url_for('static', filename=image_path))
    else:
        return render_template(template)


app.register_blueprint(bp, url_prefix='/JACKS')

if __name__ == '__main__':
    app.run()
