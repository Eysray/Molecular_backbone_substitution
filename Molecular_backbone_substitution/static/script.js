let RDKit = null;
// 初始化应用
document.addEventListener('DOMContentLoaded', function() {
    initRDKit().then(() => {
        document.getElementById('generateBtn').addEventListener('click', generateDerivatives);
        drawInputMolecule();
    });
});

async function initRDKit() {
    if (typeof window.initRDKitModule === 'function') {
        try {
            RDKit = await window.initRDKitModule();
            console.log('RDKit initialized');
        } catch (e) {
            console.error('RDKit init failed:', e);
        }
    }
}

function drawInputMolecule() {
    const smiles = document.getElementById('smilesInput').value.trim();
    drawMolecule('inputMolContainer', smiles);
}

function drawMolecule(containerId, smiles) {
    const container = document.getElementById(containerId);
    container.innerHTML = '';
    if (!smiles) return;
    // 使用RDKit
    if (RDKit) {
        try {
            const mol = RDKit.get_mol(smiles);
            if (mol) {
                container.innerHTML = mol.get_svg();
                mol.delete();
                return;
            }
        } catch (e) {
            console.error('RDKit drawing error:', e);
        }
    }
}

async function generateDerivatives() {
    drawInputMolecule();
    const btn = document.getElementById('generateBtn');
    btn.disabled = true;
    btn.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Generating...';
    
    const container = document.getElementById('derivativesContainer');
    container.innerHTML = '';
    
    const smiles = document.getElementById('smilesInput').value.trim();
    const nDerivatives = document.getElementById('derivativeCount').value.trim();

    try {
        // 准备发送的数据
        const data = { smiles: smiles, nDerivatives: nDerivatives }; 
        console.log(data);

        // 发送 POST 请求（返回 Response 对象）
        const response = await fetch('/process_data', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify(data),
        });

        // 解析响应为 JSON（需等待解析完成）
        const result = await response.json();

        // 前端处理返回数据（更新页面）
        displayDerivatives(result.derivatives);
        document.getElementById('derivativeCountBadge').textContent = `${result.derivatives.length} generated`;
        console.log(result.derivatives)
        updatePropertyCharts(result.derivatives);

    } catch (error) {
        // 统一处理错误
        console.error('请求失败:', error);
    }
    btn.disabled = false;
    btn.innerHTML = 'Generate Derivatives';
    
}

function displayDerivatives(derivatives) {
    const container = document.getElementById('derivativesContainer');
    derivatives = Array.from(derivatives); // 我添加
    derivatives.forEach((derivative, index) => {
        const div = document.createElement('div');
        div.className = 'molecule-container';
        div.innerHTML = `
            <div class="row">
                <div class="col-md-3">
                    <div id="mol-${index}" class="molecule-view"></div>
                </div>
                <div class="col-md-9">
                    <h5>Derivative ${index + 1} <span class="badge bg-info">${derivative.smiles}</span></h5>
                    <div class="mb-2">
                        <strong>SMILES:</strong> <code>${derivative.smiles}</code>
                        ${renderProperties(derivative.properties)}
                    </div>
                </div>
            </div>
        `;              
        container.appendChild(div);
        
        setTimeout(() => drawMolecule(`mol-${index}`, derivative.smiles), 100);
    });
}

function renderProperties(properties) {
    return `
        <span class="badge bg-primary property-badge">MW: ${properties.MW.toFixed(2)}</span>
        <span class="badge bg-success property-badge">LogP: ${properties.LOGP.toFixed(2)}</span>
        <span class="badge bg-danger property-badge">HBD: ${properties.HBD}</span>
        <span class="badge bg-warning text-dark property-badge">HBA: ${properties.HBA}</span>
        <span class="badge bg-info property-badge">TPSA: ${properties.TPSA.toFixed(2)}</span>
    `;
}

// 更新性质分布图表
function updatePropertyCharts(derivatives) {
    if (!derivatives || derivatives.length === 0) return;
    
    const propertyData = {
        mw: derivatives.map(d => d.properties.MW),
        logp: derivatives.map(d => d.properties.LOGP),
        hbd: derivatives.map(d => d.properties.HBD),
        hba: derivatives.map(d => d.properties.HBA)
    };
    console.log(propertyData)
    updateChart('mwChart', 'Molecular Weight Distribution', propertyData.mw, 'MW (Da)');
    updateChart('logpChart', 'LogP Distribution', propertyData.logp, 'LogP');
    updateChart('hbdChart', 'HBD Distribution', propertyData.hbd, 'HBD Count');
    updateChart('hbaChart', 'HBA Distribution', propertyData.hba, 'HBA Count');
}

// 更新单个图表
function updateChart(canvasId, title, inputdata, xLabel) {
    const ctx = document.getElementById(canvasId).getContext('2d');
    
    // 如果图表已存在，销毁它
    if (window[canvasId + 'Chart']) {
        window[canvasId + 'Chart'].destroy();
    }
    const length = inputdata.length;
    const labels = Array.from({ length }, (_, index) => index + 1)
    // 创建新图表
    window[canvasId + 'Chart'] = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: labels,
            datasets: [{
                labels: xLabel,
                data: inputdata,
                backgroundColor: 'rgba(54, 162, 235, 0.5)',
                borderColor: 'rgba(54, 162, 235, 1)',
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            plugins: {
                title: {
                    display: true,
                    text: title
                },
                legend: {
                    display: false
                }
            },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'number of derivatives'
                    }
                },
                y: {
                    title: {
                        display: true,
                        text: xLabel
                    },
                    beginAtZero: true
                }
            }
        }
    });
}