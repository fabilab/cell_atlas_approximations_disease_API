# app.py
"""
Web application supporting the cell atlas approximation disease API
"""
from flask import Flask
from flask_cors import CORS
from routes.differential_celltype_abundance import diff_celltype_abundance_bp
from routes.differential_gene_expression import diff_gene_expression_bp
from routes.get_metadata import get_metadata_bp

app = Flask(__name__)
with open('secret_key.txt') as f:
    app.config['SECRET_KEY'] = f.read()

# Cross-origin request handler
# FIXME: do not open to arbitrary requests obviously
CORS(app, resources={r"/*": {"origins": "*"}})

# Register blueprints
app.register_blueprint(diff_celltype_abundance_bp)
app.register_blueprint(diff_gene_expression_bp)
app.register_blueprint(get_metadata_bp)

if __name__ == '__main__':
    # when running locally use 127.0.0.1
    # host = '127.0.0.1'
    # when deployment, use 0.0.0.0
    host = '0.0.0.0'
    app.run(host='127.0.0.1', port=5000)
