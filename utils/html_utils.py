# HTML Utils

class HTMLUtils:
    def __init__(self, filename):
        self.filename = filename
        self.open(filename)

    def open(self, filename):
        self.fp = open(filename, 'w')
        header = '''
<html lang="en-US">
<head>
<style>
.row {
  display: flex;
}
.column {
  flex: 50%;
}
</style>
</head>
<body>
        '''
        self.fp.write('%s\n' % header)

    def add_header(self, txt):
        self.fp.write('<h1>%s</h1>\n' % txt)

    def start_div(self):
        self.fp.write('<div>\n')

    def end_div(self):
        self.fp.write('</div>\n')

    def add_text(self, txt):
        self.fp.write('%s<br>' % txt)

    def add_images(self, image_list):
        self.fp.write('<div class="row">\n')
        for image in image_list:
            self.fp.write('<div class="column">\n')
            top_txt = image.get('title', None)
            if top_txt is not None:
                self.fp.write('<center>%s</center><br>' % top_txt)
            self.fp.write('<img src="%s"><br>\n' % image.get('img', ''))
            bottom_txt = image.get('subtitle', None)
            if bottom_txt is not None:
                self.fp.write('<center>%s</center>' % bottom_txt)
            self.fp.write('</div>\n')
        self.fp.write('</div>\n')

    def add_image(self, image_name):
        self.fp.write('<img src="%s"><br>\n' % image_name)

    def close(self):
        self.fp.write('</body>\n')
        self.fp.close()
        self.fp = None

