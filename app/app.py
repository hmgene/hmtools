from flask import Flask
from flask import render_template

app = Flask(__name__)

@app.route("/")
def root(name=None):
    return render_template('sankey.html', name=name)

@app.route("/data")
def data(name=None):
    return 'data'

if __name__ == "__main__":
	app.run(debug=True)
