from flask import Flask, jsonify, render_template
import pandas as pd
import numpy as np
import sys
import json


app = Flask(__name__)
 
def prt(*a):
    # Here a is the array holding the objects
    # passed as the argument of the function
    print(*a, file=sys.stderr)
 
 

@app.route("/")
def root(name=None):
	df = pd.read_csv('static/horizon.csv');
	x=df.groupby(by="sample"),sort_values(by="length", ascending=True)["counts"].values.tolist();
	y={ "k":"k", "v":x }
	return render_template('horizon.html',data=json.dumps(y));
    #return render_template('sankey.html', name=name)

@app.route("/data")
def data(name=None):
    return 'data'

if __name__ == "__main__":
	df = pd.read_csv('static/horizon.csv');
	#g=df.groupby(by="sample"),sort_values(by="length", ascending=True)["counts"].values.tolist();
	print(df.groupby(by="sample")["counts"].apply(list))
	print(df.groupby(by="sample")["length"].apply(list))
	#y={ 'k':1, 'v':x}
	#app.run(debug=True)
