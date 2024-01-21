from flask import Flask, render_template, request, redirect, url_for,jsonify, session
from flask_sqlalchemy import SQLAlchemy 
from sqlalchemy import create_engine, MetaData, Table, Column, String, Float, Integer,or_,TEXT,text
from werkzeug.security import check_password_hash, generate_password_hash
import math
from functools import wraps
import jwt
from db import mysql_host,mysql_port,mysql_user,mysql_password,mysql_database
app = Flask(__name__)

app.config['SQLALCHEMY_DATABASE_URI'] = f'mysql+mysqlconnector://{mysql_user}:{mysql_password}@{mysql_host}/{mysql_database}'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['JWT_SECRET_KEY'] = 'your_jwt_secret_key'
app.secret_key = 'your_secret_key'

db = SQLAlchemy(app)


metadata = MetaData()

# Define the 'lingads' table
lingads = Table(
    'lingads', metadata,
    Column('pdb_code', String(255), primary_key=True),
    Column('resolution', Float),
    Column('release_year', Integer),
    Column('binding_data', String(20)),
    Column('logKdKi', String(20)),
    Column('reference', String(255)),
    Column('ligand_name', String(50)),
    extend_existing=True
)

# Define the 'proteins' table
proteins = Table(
    'proteins', metadata,
    Column('pdb_code', String(255), primary_key=True),
    Column('release_year', Integer),
    Column('uniprot_id', String(20)),
    Column('protein_name', String(255)),
    extend_existing=True
)

# Define the 'users' table
users = Table(
    'users', metadata,
    Column('id', Integer, primary_key=True),
    Column('username', String(80), unique=True, nullable=False),
    Column('password', String(255), nullable=False),
    extend_existing=True
)

lingadInfo = Table (
    'lingad_infos', metadata,
    Column('id', Integer, primary_key=True),
    Column('code', String(10), unique=True),
    Column('name', String(255)),
    Column('identifiers', String(255)),
    Column('formula', String(255)),
    Column('molecular_weight', String(255)),
    Column('lingad_type', String(255)),
    Column('additional_info', TEXT),
)

proteinInfo = Table (
    'protein_infos', metadata,
    Column('id', Integer, primary_key=True),
    Column('name', String(255)),
    Column('additional_info', TEXT),
)

# Create an engine and bind it to the metadata
engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'], echo=True)
metadata.create_all(engine)
folder_path = './upload'

class Lingads(db.Model):
    __table__ = lingads

class Proteins(db.Model):
    __table__ = proteins

class User(db.Model):
    __table__ = users

class LingadInfo(db.Model):
    __table__ = lingadInfo
    
class ProteinInfo(db.Model):
    __table__ = proteinInfo
    
def login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if 'token' not in session:
            return redirect(url_for('login'))
        try:
            jwt.decode(session['token'], app.config['JWT_SECRET_KEY'], algorithms=['HS256'])
        except jwt.ExpiredSignatureError:
            return redirect(url_for('login'))
        return f(*args, **kwargs)
    return decorated_function

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        user = User.query.filter_by(username=username).first()
        if user and check_password_hash(user.password, password):
            token = jwt.encode({'username': username}, app.config['JWT_SECRET_KEY'], algorithm='HS256')
            session['token'] = token
            return redirect(url_for('home'))
        else:
            return render_template('login.html', error='Invalid credentials')

    return render_template('login.html')

@app.route('/logout')
def logout():
    session.pop('token', None)  # Remove the token from the session
    return redirect(url_for('login'))

@app.route('/dataset', methods=['GET'])
@login_required
def dataset():
    try:
        protein_name = request.args.get('protein_name')
        code = request.args.get('code')
        page = int(request.args.get('page', 1))
        per_page = 20
        if protein_name is  None:
            protein_name = ""
            
        if code is  None:
            code = ""
        result = (
            db.session.query(lingads,proteins)
            .join(proteins, lingads.c.pdb_code == proteins.c.pdb_code)
            .filter( 
                or_(
                Proteins.protein_name.ilike(f'%{protein_name}%'),
                Lingads.ligand_name.ilike(f'%{protein_name}%'),
                ),
                Proteins.pdb_code.ilike(f'%{code}%'),
            )
            .limit(per_page).offset((page - 1) * per_page).all()
        )
        lingads_list = [
            {
                'pdb_code': row.pdb_code,
                'resolution': row.resolution,
                'logKdKi' : row.logKdKi,
                'release_year': row.release_year,
                'binding_data': row.binding_data,
                'reference': row.reference,
                'ligand_name': row.ligand_name,
                'protein_name': row.protein_name,
            }
            for row in result
        ] 
        count =(db.session.query(lingads,proteins)
            .join(proteins, lingads.c.pdb_code == proteins.c.pdb_code)
            .filter( or_(
            Proteins.protein_name.ilike(f'%{protein_name}%'),
            Lingads.ligand_name.ilike(f'%{protein_name}%')
        )).count())
     
        total_page = math.ceil(count/per_page)
        next_page = 0 
        prev_page = 0
        if total_page > page:
           next_page = page +1
           
        if page>1:
           prev_page = page -1
        return render_template('dataset.html', lingads_list=lingads_list,total=count,pname=protein_name,npage=next_page,ppage=prev_page,page=page,code=code)
    except Exception as e:
        return jsonify({'error': str(e)})

@app.route('/ligand/<string:code>', methods=['GET'])
@login_required
def ligand_info(code):
    lingadInfo = LingadInfo.query.filter_by(code=code).first()
    if lingadInfo:
        return render_template('ligand_info.html',ligand_info=lingadInfo)
    else:
        return render_template('dataset.html', error='Invalid code')

@app.route('/add/ligand', methods=['GET','POST'])
@login_required
def ligand():
    if request.method == 'POST':
        code = request.form['code']
        name = request.form['name']
        identifiers = request.form['identifiers']
        formula = request.form['formula']
        molecular_weight = request.form['molecular_weight']
        lingad_type = request.form['lingad_type']
        additional_info = request.form['additional_info']
        
        new_data = LingadInfo(
            code=code,
            name=name,
            identifiers=identifiers,
            formula=formula,
            molecular_weight=molecular_weight,
            lingad_type=lingad_type,
            additional_info=additional_info,
        )

        db.session.add(new_data)
        db.session.commit()

        return redirect(url_for('ligands'))

    return render_template('ligand_add.html')

@app.route('/ligands', methods=['GET'])
@login_required
def ligands():
    ligand_list = LingadInfo.query.all()
    return render_template('ligand_list.html', ligands=ligand_list)


from rdkit import Chem
from openmm.app import PDBFile
from pdbfixer import PDBFixer
import mdtraj as md
import nglview
from deepchem.utils.vina_utils import prepare_inputs
@app.route('/add/protein', methods=['GET','POST'])
@login_required
def protein():
    if request.method == 'POST':
        name = request.form['name']
        additional_info = request.form['additional_info']
        protein_file = request.files['protein_file']
        pocket_file = request.files['pocket_file']
        
        protein_file.save(folder_path + "/%s_protein.pdb" % (name))
        pocket_file.save(folder_path + "/%s_pocket.pdb" % (name))
        
        fixer = PDBFixer(filename=folder_path + "/%s_protein.pdb" % (name))

        PDBFile.writeFile(fixer.topology, fixer.positions, open('templates/display/%s.pdb' % (name), 'w'))
        p, m = None, None
        p, m = prepare_inputs('templates/display/%s.pdb' % (name), "C")
        if p and m:  
            Chem.rdmolfiles.MolToPDBFile(p, 'templates/display/protein_%s.pdb' % (name))
                
        protein_mdtraj = md.load_pdb( 'templates/display/protein_%s.pdb' % (name))
        p = nglview.show_mdtraj(protein_mdtraj)
        nglview.write_html("templates/display/protein_%s.html" % (name),[p])
        
        new_data = ProteinInfo(
            name=name,
            additional_info=additional_info,
        )

        db.session.add(new_data)
        db.session.commit()

        return redirect(url_for('list_protein'))

    return render_template('protein_add.html')

@app.route('/proteins', methods=['GET'])
@login_required
def list_protein():
    protein_list = ProteinInfo.query.all()
    return render_template('protein_list.html', proteins=protein_list)

import requests
import os



@app.route('/', methods=['GET','POST'])
@login_required
def home():
    lid = request.form.get('ligand')
    protein_name = request.form.get('protein')
    ligands = LingadInfo.query.all()
    
   
    if request.method == 'POST':
        selected_ligand = LingadInfo.query.filter_by(id=lid).first()
        ligand = (
            db.session.query(lingads)
            .filter( 
                or_(
                Lingads.ligand_name.ilike(f'%{selected_ligand.code}%'),
                ),
            ).first()
        )
        
        protein = (
            db.session.query(proteins)
            .filter( 
                or_(
                Proteins.name.ilike(f'%{protein_name}%'),
                ),
            ).first()
        )
        if not os.path.exists("./upload/%s_pocket.pdb" % protein_name):
            return render_template('home.html', ligands=ligands, lcode=lid, selected_ligand= selected_ligand, error = "File %s_pocket.pdb not found" %protein_name, protein_name=protein_name)
        if ligand is None or protein is None:
            return render_template('home.html', ligands=ligands, lcode=lid, selected_ligand= selected_ligand, error = "Not fond data in DB",protein_name=protein_name)
        url = "http://127.0.0.1:5001?ligand=%s&protein=%s" % (ligand[0],protein_name)

      
        response = requests.request("GET", url, headers={
        'Content-Type': 'application/json'
        })
        if response.status_code == 200:
            json_data = response.json()
            result_value = round(json_data['result'][0], 2)
        else:
            return render_template('home.html', ligands=ligands, lcode=lid, selected_ligand= selected_ligand, error = "Model run fail")
        return render_template('home.html', ligands=ligands, lcode=lid, selected_ligand= selected_ligand, lpdb_code = ligand[0], protein_name = protein_name,result =result_value,protein=protein)
    
    return render_template('home.html', ligands=ligands, lcode=lid)

@app.route('/display/<string:code>', methods=['GET'])
def display(code):
    return render_template('display/'+ code +'.html')

if __name__ == "__main__":
    app.run(host="0.0.0.0",port=5005)