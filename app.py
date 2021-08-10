import os
from flask import Flask, render_template, request, redirect, url_for, abort
from werkzeug.utils import secure_filename
import needleman
import gibbsSampling

app = Flask(__name__)
app.config["DEBUG"] = True

app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024  # Upload edilebilecek maksimum dosya boyutunu 1MB olarak ayarladım.
app.config['UPLOAD_EXTENSIONS'] = ['.txt']  # Sadece .txt dosyası upload edilmesini sağlıyor
app.config['UPLOAD_PATH'] = 'uploads'  # Yüklenen dosyaların "uploads" isimli dosyaya kaydedilmesini sağlıyor


@app.route("/")
def main():
    return render_template('index.html')


@app.route('/', methods=['GET', 'POST'])
def upload_files():
    uploaded_file = request.files['file']
    filename = secure_filename(uploaded_file.filename)
    if filename != '':
        file_ext = os.path.splitext(filename)[1]
        if file_ext not in app.config['UPLOAD_EXTENSIONS']:
            abort(400)
        uploaded_file.save(os.path.join(app.config['UPLOAD_PATH'],
                                        "dnaSequence1.txt"))  # Yüklenen dosyanın ismini istediğimiz gibi değiştiriyor

    if request.method == "POST":
        # getting input with name = fname in HTML form
        motifL = request.form.get("motifLength")

        needleman.needle_calculations()  # Hesaplama yapılacak fonksiyona yönlendirme   #Burada kaldık, yüklenen text dosyasını tek bir isim olacak şekilde saklayacak. İşlem ve sonuç için başka fonksiyona aktaracak.

        gibbsSampling.goToGibbs(motifL)

    return content()


@app.route("/")
def content():
    with open('outputs/seqAlignmentOutput.txt', 'r') as f, open('outputs/motifFindingOutput.txt', 'r') as f2:
        return render_template('after_upload.html', seqAlignmentFinal=f.read(), motifFindingOutput=f2.read())


if __name__ == "__main__":
    app.run()(Debug=true)
