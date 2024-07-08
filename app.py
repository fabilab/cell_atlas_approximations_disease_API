# app.py

from flask import Flask
from flask_cors import CORS
from routes.differential_celltype_abundance import diff_celltype_abundance_bp
from routes.differential_gene_expression import diff_gene_expression_bp

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "http://localhost:3000"}})

# Register blueprints
app.register_blueprint(diff_celltype_abundance_bp)
app.register_blueprint(diff_gene_expression_bp)

if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=False)