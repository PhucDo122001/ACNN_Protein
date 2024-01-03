import deepchem as dc
import os

import numpy as np
# import tensorflow as tf

import matplotlib.pyplot as plt

from deepchem.molnet import load_pdbbind
from deepchem.models import AtomicConvModel
from deepchem.feat import AtomicConvFeaturizer

f1_num_atoms = 100  # maximum number of atoms to consider in the ligand
f2_num_atoms = 1000  # maximum number of atoms to consider in the protein
max_num_neighbors = 12  # maximum number of spatial neighbors for an atom

acf = AtomicConvFeaturizer(frag1_num_atoms=f1_num_atoms,
                      frag2_num_atoms=f2_num_atoms,
                      complex_num_atoms=f1_num_atoms+f2_num_atoms,
                      max_num_neighbors=max_num_neighbors,
                      neighbor_cutoff=4)
tasks, datasets, transformers = load_pdbbind(featurizer=acf,
                                             save_dir='.',
                                             data_dir='.',
                                             pocket=True,
                                             reload=False,
                                             set_name='core')
class MyTransformer(dc.trans.Transformer):
  def transform_array(x, y, w, ids):
    kept_rows = x != None
    X_reshaped = np.squeeze(x[kept_rows], axis=0)  # Remove the first dimension
    Y_reshaped = np.squeeze(y[kept_rows], axis=0)  # Remove the first dimension
    W_reshaped = np.squeeze(w[kept_rows], axis=0)  # Remove the first dimension
    ID_reshaped = np.squeeze(ids[kept_rows], axis=0)  # Remove the first dimension
    return X_reshaped, Y_reshaped, W_reshaped, ID_reshaped

datasets = [d.transform(MyTransformer) for d in datasets]
train, val, test = datasets
acm = AtomicConvModel(n_tasks=1,
                      frag1_num_atoms=f1_num_atoms,
                      frag2_num_atoms=f2_num_atoms,
                      complex_num_atoms=f1_num_atoms+f2_num_atoms,
                      max_num_neighbors=max_num_neighbors,
                      batch_size=4,
                      layer_sizes=[32, 32, 16],
                      learning_rate=0.003,
                      )
losses, val_losses = [], []
#max_epochs = 50
max_epochs = 1 

metric = dc.metrics.Metric(dc.metrics.score_function.rms_score)
step_cutoff = len(train)//12
def val_cb(model, step):
    if step%step_cutoff!=0:
        return
    val_losses.append(model.evaluate(val, metrics=[metric])['rms_score']**2)  # L2 Loss
    losses.append(model.evaluate(train, metrics=[metric])['rms_score']**2)  # L2 Loss
    f, ax = plt.subplots()
    ax.scatter(range(len(losses)), losses, label='train loss')
    ax.scatter(range(len(val_losses)), val_losses, label='val loss')
    plt.legend(loc='upper right')

acm.fit(train, nb_epoch=max_epochs, max_checkpoints_to_keep=1,
                callbacks=[val_cb])
import pandas as pd
index_labels_file = "./v2013-core/pdbbind_v2013_core.csv"
data_folder = "./v2013-core/"
df = pd.read_csv(index_labels_file)
pdbs = df.pdb_id.tolist()
protein_files = [
                os.path.join(data_folder, pdb, "%s_pocket.pdb" % pdb)
                for pdb in pdbs[:2]
            ]
     
ligand_files = [
            os.path.join(data_folder, pdb, "%s_ligand.sdf" % pdb)
            for pdb in pdbs[:2]
        ]

from flask import Flask, flash, redirect, url_for,request,jsonify

app = Flask(__name__)
app.secret_key = "your_secret_key"  # Required for flash messages

@app.route('/', methods=['GET'])
def login():
    try:
        ligand_pdbcode = request.args.get("ligand")
        protein_pdbcode = request.args.get("protein")
        if ligand_pdbcode is None or protein_pdbcode is None:
            ligand_pdbcode ="3zsx"    
            protein_pdbcode ="3zsx"    
        protein_files[0] = os.path.join(data_folder, protein_pdbcode, "%s_pocket.pdb" % protein_pdbcode)
        ligand_files[0] =os.path.join(data_folder, ligand_pdbcode, "%s_ligand.sdf" % ligand_pdbcode) 
        f2 = acf.featurize(list(zip(ligand_files, protein_files)))
        input = dc.data.DiskDataset.from_numpy(f2)
        npa = acm.predict(input)
        response_data = {
            "result": npa[0][0].tolist()
        }

        return jsonify(response_data)
    except Exception as e:
        print(e)
        return jsonify({"error": "error"}), 400

if __name__ == "__main__":
    app.run(host="0.0.0.0",port=5001)