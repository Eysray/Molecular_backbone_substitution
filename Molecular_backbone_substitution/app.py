from flask import Flask, render_template, request, jsonify
from toolkit import replacecore_many as rm
from toolkit import add_smi_prop as addp
import os

app = Flask(__name__)

@app.route('/process_data', methods=['POST'])
def process_data():
        # 接收前端发送的数据
    data = request.get_json()
    smiles = data.get('smiles')# 假设前端发送JSON数据
    print(smiles)
    nDerivatives = data.get('nDerivatives')
    nDerivatives = int(nDerivatives)
    print(nDerivatives)
    dataset_path = r'data\merged.csv'
    # 模拟数据处理（例如反转字符串）
    derivatives = rm(smiles, dataset_path, nDerivatives)
    derivatives = list(set(derivatives)) #去重,乱序
    output = []
    for derivative in derivatives:
        output.append(addp(derivative))
    '''
    output结构:
    [
        {
            'smiles': smi: str, 
            'properties':
                {
                    'prop_name': prop_value, 
                },
        },
    ]
    '''
    # 返回处理结果给前端
    return jsonify({'derivatives': output})


@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=False)