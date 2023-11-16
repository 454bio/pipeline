import sys
import html_utils

in_filename = 'run.log'
out_filename = 'run.html'

argcc = 1
argc = len(sys.argv)
while argcc < argc:
    if sys.argv[argcc] == '-i':
        argcc += 1
        in_filename = sys.argv[argcc]
    if sys.argv[argcc] == '-o':
        argcc += 1
        out_filename = sys.argv[argcc]
    argcc += 1

html = html_utils.HTMLUtils(out_filename)
html.add_header('Analysis Report')

html.start_div()
html.add_text('<h2>phase parameter estimation</h2>')
html.add_images([{'img':'err.png', 'subtitle': 'phase fit error'}])
html.add_images([
    {'img':'ie.png', 'subtitle':'incomplete extension estimates'},
    {'img':'cf.png', 'subtitle':'carry-forward estimates'},
    {'img':'dr.png', 'subtitle':'signal loss estimates'}
])
html.end_div()

html.start_div()
html.add_text('<h2>Alignment</h2>')
html.add_images([
    {'img':'coverage.png', 'subtitle':'coverage from all reads'},
    {'img':'filtered_coverage.png', 'subtitle':'filtered read coverage'},
    {'img':'coverage_top_50Q10.png', 'subtitle':'coverage from top 50 Q10 reads'}
])
html.end_div()

html.start_div()
html.add_text('<h2>Run Log</h2>')
html.add_text('<code><pre>')
with open(in_filename) as f:
    lines = f.readlines()
    for line in lines:
        html.add_text(line.rstrip())
html.add_text('</pre></code>')
html.end_div()

html.close()

