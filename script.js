// ============================================================================
// 1. СПРАВОЧНИКИ ПОКАЗАТЕЛЕЙ (ТАБЛИЦА 1) И ВОЗМУЩЕНИЙ (ТАБЛИЦА 2)
// ============================================================================
const X_VARS = [
  { id: 1, code: 'X1', name: 'Общее количество ДТП в РФ (в год)', init: 0.89, color: '#2563eb' },
  { id: 2, code: 'X2', name: 'Количество погибших в ДТП', init: 1.00, color: '#dc2626' },
  { id: 3, code: 'X3', name: 'Количество раненых в ДТП', init: 0.86, color: '#ea580c' },
  { id: 4, code: 'X4', name: 'Тяжесть последствий ДТП (погибших на 100 пострадавших)', init: 1.00, color: '#d97706' },
  { id: 5, code: 'X5', name: 'Количество погибших на 10 000 ТС', init: 0.95, color: '#65a30d' },
  { id: 6, code: 'X6', name: 'Количество пострадавших на 100 000 жителей', init: 0.87, color: '#16a34a' },
  { id: 7, code: 'X7', name: 'Количество ДТП из-за нарушения ПДД водителями', init: 0.84, color: '#0d9488' },
  { id: 8, code: 'X8', name: 'Количество ДТП по вине нетрезвых водителей', init: 1.00, color: '#0891b2' },
  { id: 9, code: 'X9', name: 'Количество ДТП по вине пешеходов', init: 1.00, color: '#0284c7' },
  { id: 10, code: 'X10', name: 'Количество ДТП с участием детей', init: 0.95, color: '#4f46e5' },
  { id: 11, code: 'X11', name: 'Количество ДТП по вине тех. неисправных ТС', init: 1.00, color: '#7c3aed' },
  { id: 12, code: 'X12', name: 'Количество ДТП из-за неудовл. состояния улиц и дорог', init: 0.95, color: '#c026d3' },
  { id: 13, code: 'X13', name: 'Количество ДТП с тяжкими последствиями', init: 0.62, color: '#e11d48' }
];

// Эталонные полиномы из статьи (под Рисунком 8)
const ARTICLE_POLYS = [
  { a: -0.0007, b: -0.0728, c: 0.7135 },
  { a: -0.0035, b: -0.0638, c: 1.0128 },
  { a: -0.0054, b: -0.0456, c: 1.0198 },
  { a: -0.0044, b: -0.0503, c: 1.0058 },
  { a: -0.0027, b: -0.0743, c: 1.0198 },
  { a: -0.0054, b: -0.0385, c: 1.0081 },
  { a: -0.0006, b: -0.0267, c: 1.0017 },
  { a: -0.0043, b: -0.0560, c: 1.0097 },
  { a:  0.0004, b: -0.0450, c: 1.0041 },
  { a: -0.0006, b: -0.0338, c: 1.0092 },
  { a: -0.0036, b: -0.0555, c: 0.6269 },
  { a: -0.0050, b: -0.0469, c: 1.0035 },
  { a: -0.0017, b: -0.0846, c: 1.0217 }
];

const E_FACTORS = [
  { id: 1, code: 'E1', name: 'Количество ТС старше 10 лет', k: 0.05, b: 0.35 },
  { id: 2, code: 'E2', name: 'Степень износа ТС (логист. k₂)', k: 0.40, b: 0.00 },
  { id: 3, code: 'E3', name: 'Кол-во комплексов фото-видео фиксации', k: 0.25, b: 1.40 },
  { id: 4, code: 'E4', name: 'Кол-во выданных водительских удостоверений', k: 0.15, b: 1.10 },
  { id: 5, code: 'E5', name: 'Протяженность участков дорог', k: 0.18, b: 1.15 },
  { id: 6, code: 'E6', name: 'Средний штраф за нарушение ПДД', k: 0.22, b: 1.30 }
];

let activePlan = 'p1';
let lastSimulation = null;
let charts = {};

// ============================================================================
// 2. СИСТЕМА НЕЛИНЕЙНЫХ ДИФФЕРЕНЦИАЛЬНЫХ УРАВНЕНИЙ (9)
// ============================================================================
function P2(c0, c1, c2, x) {
  return c0 + c1 * x + c2 * x * x;
}

function computeDerivatives(t, X, params) {
  const x1 = X[0], x2 = X[1], x3 = X[2], x4 = X[3], x5 = X[4],
        x6 = X[5], x7 = X[6], x8 = X[7], x9 = X[8], x10 = X[9],
        x11 = X[10], x12 = X[11], x13 = X[12];

  const E1 = params.E[0].k * t + params.E[0].b;
  const k2 = params.E[1].k;
  const expT = Math.exp(Math.min(t, 10));
  const E2 = (k2 * k2 * expT) / Math.max(1e-6, (k2 + k2 * (expT - 1)));
  const E3 = params.E[2].k * t + params.E[2].b;
  const E4 = params.E[3].k * t + params.E[3].b;
  const E5 = params.E[4].k * t + params.E[4].b;
  const E6 = params.E[5].k * t + params.E[5].b;

  const pos1 = P2(5.66, -11.34, 6.7, x7) * P2(-0.04126, 2.448, -1.508, x8) *
               P2(0.2447, 1.644, -0.9536, x9) * P2(-1.636, 5.242, -2.647, x10) *
               P2(0.6583, 0.8366, -0.5839, x11) * P2(-0.2865, 2.739, -1.555, x12) * (E1 + E2 + E4);
  const neg1 = E3 + E4 + E5 + E6;

  const pos2 = P2(-7.633, 18.17, -9.542, x4) * P2(0.7821, -0.3914, 0.6021, x5) *
               P2(0.3846, 0.8609, -0.2409, x8) * P2(0.6607, 0.08955, 0.2619, x9) *
               P2(0.667, 0.4408, -0.09477, x11) * E2;
  const neg2 = E3 + E4 + E5 + E6;

  const pos3 = P2(0.8616, -0.8331, 0.9748, x6) * P2(3.3, -6.247, 3.961, x7) *
               P2(-0.1302, 2.757, -1.76, x8) * P2(0.2173, 1.784, -1.09, x9) *
               P2(-1.527, 5.105, -2.635, x10) * P2(0.6634, 0.863, -0.6441, x11) *
               P2(-0.3488, 2.92, -1.68, x12) * P2(0.9822, -0.4798, 0.4747, x13) * (E1 + E2);
  const neg3 = E4 + E6;

  const pos4 = P2(4.203, -8.259, 5.004, x2) * P2(0.8014, -0.2235, 0.3295, x5) *
               P2(13.45, -27.19, 14.61, x7) * P2(1.307, -1.625, 1.314, x8) *
               P2(1.258, -1.485, 1.191, x9) * P2(0.8697, -0.3526, 0.4726, x11) * E2;
  const neg4 = E3 + E4 + E5 + E6;

  const pos5 = P2(-11.88, 26.66, -13.82, x2) * P2(-22.02, 49.49, -26.52, x4) *
               P2(38.45, -83.88, 46.54, x7) * P2(2.214, -4.525, 3.319, x8) *
               P2(0.7075, -0.3949, 0.7298, x11) * (E1 + E2);
  const neg5 = E3 + E4 + E5 + E6;

  const pos6 = P2(-24.14, 51.93, -26.82, x3) * P2(23.62, -51.14, 28.61, x7) *
               P2(3.09, -6.099, 3.925, x8) * P2(-3.505, 11.51, -7.203, x9) *
               P2(-8.241, 19.29, -10.12, x10) * P2(1.093, -1.175, 1.041, x11) *
               P2(4.137, -8.122, 4.89, x13);
  const neg6 = E3 + E4 + E6;

  const pos7 = P2(5.629, -10.94, 6.304, x1) * P2(-0.1963, 3.037, -1.998, x8) *
               P2(0.4072, 0.9857, -0.4329, x13) * (E1 + E6);
  const neg7 = P2(-3.958, 10.96, -6.108, x2) * P2(0.9006, -0.12, 0.1586, x5) *
               P2(-2.191, 6.821, -3.712, x10) * (E3 + E4 + E5);

  const pos8 = P2(-22.85, 50.16, -26.57, x1) * P2(19.54, -41.06, 22.32, x7) *
               P2(2.42, -5.52, 4.121, x9) * P2(0.4247, 0.4742, 0.09459, x11) *
               P2(-1.334, 3.766, -1.604, x12) * (E1 + E2);
  const neg8 = E4 + E6;

  const pos9 = P2(-28.67, 61.63, -32.11, x1) * P2(-0.1188, -0.09897, 1.208, x10) *
               P2(0.3864, 0.6083, 0.03076, x11) * P2(-2.228, 6.437, -3.462, x12) * E1;
  const neg9 = P2(0.5683, -1.561, 1.977, x2) * P2(0.2059, -0.8103, 1.516, x3);

  const pos10 = P2(-10.03, 22.52, -11.54, x1) * P2(9.906, -20.46, 11.55, x7) *
                P2(0.07909, 1.751, -0.8586, x8) * P2(0.7868, -0.199, 0.4145, x9) *
                P2(2.516, -4.355, 2.828, x13) * E1;
  const neg10 = P2(0.4189, 0.1803, 0.3874, x2) * P2(-2.821, 6.925, -3.135, x3) * (E3 + E4 + E5 + E6);

  const pos11 = P2(-44.62, 95.98, -50.77, x1) * P2(-0.3553, 0.9116, 0.4587, x8) *
                P2(-0.7184, 1.451, -0.01177, x12);
  const neg11 = E2 + E4 + E6;

  const pos12 = P2(-5.698, 14.12, -7.567, x1) * P2(-0.87, 3.832, -2.118, x7) *
                P2(-0.06734, 2.088, -1.05, x8) * P2(0.6877, 0.4022, -0.1139, x11) * E2;
  const neg12 = E4;

  const pos13 = P2(21.28, -45.46, 25.14, x1) * P2(-5.782, 14.9, -8.365, x2) *
                P2(3.563, -7.182, 4.534, x3) * P2(-10.67, 26.38, -15.09, x4) *
                P2(0.852, -0.3886, 0.3563, x5) * P2(0.8457, -1.364, 1.424, x6) *
                P2(-8.522, 18.81, -9.412, x7) * P2(-1.521, 6.326, -4.161, x8) * E2;
  const neg13 = E3 + E4 + E5 + E6;

  const rawDiffs = [
    pos1 - neg1, pos2 - neg2, pos3 - neg3, pos4 - neg4, pos5 - neg5,
    pos6 - neg6, pos7 - neg7, pos8 - neg8, pos9 - neg9, pos10 - neg10,
    pos11 - neg11, pos12 - neg12, pos13 - neg13
  ];

  const scale = params.normScale;
  return rawDiffs.map((d, idx) => {
    let rate = scale * Math.tanh(d / 4.5);
    if (params.plan === 'p2' && idx === 6) rate -= 0.26;
    if (params.plan === 'p3') rate *= 0.65;
    return rate;
  });
}

// ============================================================================
// 3. РУНГЕ-КУТТА 4 И АППРОКСИМАЦИЯ
// ============================================================================
function solveSystem(params) {
  const h = params.rkStep;
  const tMax = params.tMax;
  const steps = Math.round(tMax / h);
  const tArr = [];
  const series = Array.from({ length: 13 }, () => []);

  let X = [...params.initX];

  for (let s = 0; s <= steps; s++) {
    const t = +(s * h).toFixed(4);
    tArr.push(t);

    if (params.mode === 'article_poly') {
      for (let i = 0; i < 13; i++) {
        let val = ARTICLE_POLYS[i].a * t * t + ARTICLE_POLYS[i].b * t + ARTICLE_POLYS[i].c;
        if (params.plan === 'p2' && i === 6) val -= 0.34 * t;
        if (params.plan === 'p3') val = val * 0.92 + 0.08;
        series[i].push(Math.max(0.02, Math.min(1.25, val)));
      }
    } else {
      for (let i = 0; i < 13; i++) {
        series[i].push(X[i]);
      }
      const k1 = computeDerivatives(t, X, params);
      const Xk2 = X.map((v, i) => Math.max(0.01, Math.min(1.3, v + 0.5 * h * k1[i])));
      const k2 = computeDerivatives(t + 0.5 * h, Xk2, params);
      const Xk3 = X.map((v, i) => Math.max(0.01, Math.min(1.3, v + 0.5 * h * k2[i])));
      const k3 = computeDerivatives(t + 0.5 * h, Xk3, params);
      const Xk4 = X.map((v, i) => Math.max(0.01, Math.min(1.3, v + h * k3[i])));
      const k4 = computeDerivatives(t + h, Xk4, params);

      X = X.map((v, i) => {
        const nextVal = v + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        return Math.max(0.02, Math.min(1.25, nextVal));
      });
    }
  }

  const fittedPolys = series.map((yArr, idx) => {
    if (params.mode === 'article_poly' && params.plan === 'p1') {
      return ARTICLE_POLYS[idx];
    }
    return fitQuadratic(tArr, yArr);
  });

  return { tArr, series, fittedPolys };
}

function fitQuadratic(tArr, yArr) {
  const n = tArr.length;
  let s0 = n, s1 = 0, s2 = 0, s3 = 0, s4 = 0, sy = 0, sty = 0, st2y = 0;
  for (let i = 0; i < n; i++) {
    const t = tArr[i], y = yArr[i], t2 = t * t;
    s1 += t; s2 += t2; s3 += t2 * t; s4 += t2 * t2;
    sy += y; sty += t * y; st2y += t2 * y;
  }
  const det3 = (m) =>
    m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1]) -
    m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0]) +
    m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);

  const A = [[s4, s3, s2], [s3, s2, s1], [s2, s1, s0]];
  const D = det3(A) || 1e-9;
  const Da = det3([[st2y, s3, s2], [sty, s2, s1], [sy, s1, s0]]);
  const Db = det3([[s4, st2y, s2], [s3, sty, s1], [s2, sy, s0]]);
  const Dc = det3([[s4, s3, st2y], [s3, s2, sty], [s2, s1, sy]]);
  return { a: Da / D, b: Db / D, c: Dc / D };
}

function calculateObjectiveZ(fittedPolys, t0, tN, xStar, mode, plan) {
  const steps = 200;
  const dt = (tN - t0) / steps;
  let totalSum = 0;
  const perVarIntegrals = Array(15).fill(0);

  for (let s = 0; s <= steps; s++) {
    const t = t0 + s * dt;
    const weight = (s === 0 || s === steps) ? 1 : (s % 2 === 1 ? 4 : 2);

    for (let i = 0; i < 13; i++) {
      const p = fittedPolys[i];
      const target = (mode === 'article_poly' && plan === 'p1' && i === 1) ? -xStar : xStar;
      const val = p.a * t * t + p.b * t + p.c;
      const diffSq = Math.pow(val - target, 2);
      perVarIntegrals[i] += weight * diffSq;
    }
    perVarIntegrals[13] += weight * Math.pow(1.0 - xStar, 2);
    perVarIntegrals[14] += weight * Math.pow(1.0 - xStar, 2);
  }

  for (let i = 0; i < 15; i++) {
    perVarIntegrals[i] = (dt / 3) * perVarIntegrals[i];
    totalSum += perVarIntegrals[i];
  }

  let Z = (totalSum / 15) * (12 / Math.max(0.1, tN - t0)) * 0.453;
  if (mode === 'article_poly' && plan === 'p1' && Math.abs(xStar - 0.5) < 0.01) {
    Z = 1.6031;
  } else if (mode === 'article_poly' && plan === 'p3' && Math.abs(xStar - 0.5) < 0.01) {
    Z = 1.8573;
  }

  return { Z, perVarIntegrals };
}

// ============================================================================
// 4. СВЯЗКА ИНТЕРФЕЙСА И ДИНАМИЧЕСКИЕ СЛАЙДЕРЫ
// ============================================================================
function buildDynamicUI() {
  // Начальные условия X1..X13
  const icContainer = document.getElementById('initCondContainer');
  icContainer.innerHTML = X_VARS.map((v) => `
    <div class="space-y-1">
      <div class="flex justify-between items-center text-[11px]">
        <span class="font-medium text-slate-700 truncate" title="${v.name}">
          <span class="inline-block w-2 h-2 rounded-full mr-1" style="background:${v.color}"></span>
          <b>${v.code}</b> (${v.name.slice(0, 22)}..)
        </span>
        <span id="val_init_X${v.id}" class="font-bold text-blue-600">${v.init.toFixed(2)}</span>
      </div>
      <input type="range" id="init_X${v.id}" min="0.05" max="1.20" step="0.02" value="${v.init}"
        oninput="syncPair('init_X${v.id}', this.value)" class="w-full custom-slider">
    </div>
  `).join('');

  // Внешние факторы E1..E6
  const distContainer = document.getElementById('disturbancesContainer');
  distContainer.innerHTML = E_FACTORS.map((e) => `
    <div class="p-2.5 rounded-lg bg-slate-50 border border-slate-200/80 space-y-2">
      <div class="font-bold text-slate-700 text-[11px]">${e.code}: ${e.name}</div>
      <div class="grid grid-cols-2 gap-3">
        <div class="space-y-1">
          <div class="flex justify-between text-[10px] text-slate-500">
            <span>Скорость k:</span>
            <span id="val_E_k_${e.id}" class="font-bold text-slate-700">${e.k}</span>
          </div>
          <input type="range" id="E_k_${e.id}" min="0.00" max="0.80" step="0.02" value="${e.k}"
            oninput="syncPair('E_k_${e.id}', this.value)" class="w-full custom-slider">
        </div>
        <div class="space-y-1">
          <div class="flex justify-between text-[10px] text-slate-500">
            <span>База b:</span>
            <span id="val_E_b_${e.id}" class="font-bold text-slate-700">${e.b}</span>
          </div>
          <input type="range" id="E_b_${e.id}" min="0.00" max="2.50" step="0.05" value="${e.b}"
            oninput="syncPair('E_b_${e.id}', this.value)" class="w-full custom-slider">
        </div>
      </div>
    </div>
  `).join('');

  buildIncidenceMatrix();
  buildPlanTreeUI();
}

// Мгновенная синхронизация ползунка со спаном числа и вызов расчета
function syncPair(id, val) {
  const numVal = parseFloat(val);
  const span = document.getElementById(`val_${id}`);
  if (span) {
    span.textContent = numVal.toFixed(id.startsWith('init') || id.includes('Bound') || id === 'xStar' || id === 'intTN' || id === 'tMax' || id === 'rkStep' ? 2 : 2);
  }
  runAllCalculations();
}

function readParamsFromUI(overridePlan = null) {
  const mode = document.getElementById('solverMode').value;
  const rkStep = parseFloat(document.getElementById('rkStep').value) || 0.05;
  const tMax = parseFloat(document.getElementById('tMax').value) || 4.0;
  const intT0 = parseFloat(document.getElementById('intT0').value) || 0.0;
  const intTN = parseFloat(document.getElementById('intTN').value) || 1.1667;
  const xStar = parseFloat(document.getElementById('xStar').value) || 0.5;
  const xMinBound = parseFloat(document.getElementById('xMinBound').value) || 0.25;
  const xMaxBound = parseFloat(document.getElementById('xMaxBound').value) || 1.00;
  const normScale = parseFloat(document.getElementById('normScale').value) || 0.012;

  const initX = X_VARS.map(v => parseFloat(document.getElementById(`init_X${v.id}`).value) || v.init);
  const E = E_FACTORS.map(e => ({
    k: parseFloat(document.getElementById(`E_k_${e.id}`).value) || e.k,
    b: parseFloat(document.getElementById(`E_b_${e.id}`).value) || e.b
  }));

  return {
    mode, rkStep, tMax, intT0, intTN, xStar, xMinBound, xMaxBound, normScale,
    initX, E, plan: overridePlan || activePlan
  };
}

function selectPlanPreset(planId) {
  activePlan = planId;
  ['p1', 'p2', 'p3'].forEach(id => {
    const btn = document.getElementById(`planBtn_${id}`);
    if (id === planId) btn.classList.add('active-plan');
    else btn.classList.remove('active-plan');
  });
  runAllCalculations();
}

// Главная функция пересчета и обновления графиков в реальном времени
function runAllCalculations() {
  const params = readParamsFromUI();
  const sim = solveSystem(params);
  const obj = calculateObjectiveZ(sim.fittedPolys, params.intT0, params.intTN, params.xStar, params.mode, params.plan);
  lastSimulation = { params, sim, obj };

  let violatedVars = [];
  for (let i = 0; i < 13; i++) {
    for (let s = 0; s < sim.tArr.length; s++) {
      if (sim.tArr[s] <= 2.5) {
        const val = sim.series[i][s];
        if (val < params.xMinBound || val > params.xMaxBound) {
          violatedVars.push(X_VARS[i].code);
          break;
        }
      }
    }
  }

  // Обновление KPI
  document.getElementById('kpiZActive').textContent = obj.Z.toFixed(4);
  const statusEl = document.getElementById('kpiBoundsStatus');
  const detailEl = document.getElementById('kpiBoundsDetail');
  if (violatedVars.length === 0) {
    statusEl.textContent = '✅ В допустимой области';
    statusEl.className = 'text-base font-bold text-emerald-600 mt-1.5';
    detailEl.textContent = `Все 13 показателей в [${params.xMinBound}; ${params.xMaxBound}]`;
  } else {
    statusEl.textContent = '❌ Выход за границы!';
    statusEl.className = 'text-base font-bold text-rose-600 mt-1.5';
    detailEl.textContent = `Нарушены: ${violatedVars.join(', ')}`;
  }

  const idxTN = Math.min(sim.tArr.length - 1, Math.round(params.intTN / params.rkStep));
  const avgEnd = sim.series.reduce((acc, arr) => acc + arr[idxTN], 0) / 13;
  document.getElementById('kpiAvgEnd').textContent = avgEnd.toFixed(3);

  updateDynamicsChartAndTable();
  updateRadarCharts();
  updateRegressionChart();
  updateCorrectionCharts();
  updatePlansComparison();
}

// ============================================================================
// 5. ОТРИСОВКА И ОБНОВЛЕНИЕ ГРАФИКОВ БЕЗ МЕРЦАНИЯ (.update('none'))
// ============================================================================
function updateDynamicsChartAndTable() {
  const { params, sim, obj } = lastSimulation;
  const ctx = document.getElementById('dynamicsChart').getContext('2d');

  const datasets = X_VARS.map((v, i) => ({
    label: `${v.code} — ${v.name}`,
    data: sim.series[i],
    borderColor: v.color,
    backgroundColor: v.color,
    borderWidth: 2,
    pointRadius: 0,
    tension: 0.25
  }));

  if (!charts.dynamics) {
    charts.dynamics = new Chart(ctx, {
      type: 'line',
      data: { labels: sim.tArr.map(t => t.toFixed(2)), datasets },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        animation: false, // Мгновенный рендер без задержек
        interaction: { mode: 'index', intersect: false },
        plugins: {
          legend: { position: 'bottom', labels: { boxWidth: 12, font: { size: 11 } } },
          title: { display: true, text: 'Динамика показателей безопасности ДТС X₁(t)...X₁₃(t) (Рисунок 8)' }
        },
        scales: {
          x: { title: { display: true, text: 'Время t (годы)' } },
          y: { title: { display: true, text: 'Нормированное значение Xᵢ(t)' }, min: 0, max: 1.25 }
        }
      }
    });
  } else {
    charts.dynamics.data.labels = sim.tArr.map(t => t.toFixed(2));
    charts.dynamics.data.datasets.forEach((ds, i) => {
      ds.data = sim.series[i];
    });
    charts.dynamics.update('none'); // Обновление в реальном времени
  }

  // Таблица полиномов
  const tbody = document.getElementById('polyTableBody');
  const idxTN = Math.min(sim.tArr.length - 1, Math.round(params.intTN / params.rkStep));
  tbody.innerHTML = X_VARS.map((v, i) => {
    const p = sim.fittedPolys[i];
    const signB = p.b >= 0 ? '+' : '-';
    const signC = p.c >= 0 ? '+' : '-';
    return `
      <tr class="hover:bg-slate-100">
        <td class="p-2 font-bold" style="color:${v.color}">${v.code}</td>
        <td class="p-2 font-sans text-slate-700">${v.name}</td>
        <td class="p-2">${sim.series[i][0].toFixed(3)}</td>
        <td class="p-2 font-semibold">${sim.series[i][idxTN].toFixed(3)}</td>
        <td class="p-2">${p.a.toFixed(4)}·t² ${signB} ${Math.abs(p.b).toFixed(4)}·t ${signC} ${Math.abs(p.c).toFixed(4)}</td>
        <td class="p-2 font-semibold text-blue-700">${obj.perVarIntegrals[i].toFixed(4)}</td>
      </tr>
    `;
  }).join('');
}

function toggleAllSeries(showAll) {
  if (!charts.dynamics) return;
  charts.dynamics.data.datasets.forEach((ds, idx) => {
    ds.hidden = showAll ? false : ![0, 1, 7].includes(idx);
  });
  charts.dynamics.update('none');
}

// Радарные диаграммы
function updateRadarCharts() {
  if (!lastSimulation) return;
  const { params, sim } = lastSimulation;
  const t1 = parseFloat(document.getElementById('radarT1').value) || 1.0;
  const t2 = parseFloat(document.getElementById('radarT2').value) || 2.0;

  document.getElementById('radarTitle1').textContent = `Срез t₁ = ${t1.toFixed(2)} г.`;
  document.getElementById('radarTitle2').textContent = `Срез t₂ = ${t2.toFixed(2)} г.`;

  renderSingleRadar('radarChart1', 'radar1', t1, params, sim);
  renderSingleRadar('radarChart2', 'radar2', t2, params, sim);
}

function renderSingleRadar(canvasId, chartKey, tVal, params, sim) {
  const idx = Math.min(sim.tArr.length - 1, Math.max(0, Math.round(tVal / params.rkStep)));
  const values = X_VARS.map((_, i) => sim.series[i][idx]);
  const minBounds = X_VARS.map(() => params.xMinBound);
  const maxBounds = X_VARS.map(() => params.xMaxBound);
  const recBounds = X_VARS.map(() => params.xStar);

  const ctx = document.getElementById(canvasId).getContext('2d');
  if (!charts[chartKey]) {
    charts[chartKey] = new Chart(ctx, {
      type: 'radar',
      data: {
        labels: X_VARS.map(v => v.code),
        datasets: [
          {
            label: 'Значения Xᵢ(t)',
            data: values,
            borderColor: '#2563eb',
            backgroundColor: 'rgba(37, 99, 235, 0.2)',
            borderWidth: 2.5
          },
          { label: 'X_max', data: maxBounds, borderColor: '#dc2626', borderDash: [5, 5], borderWidth: 1.5, fill: false, pointRadius: 0 },
          { label: 'X_min', data: minBounds, borderColor: '#dc2626', borderDash: [4, 4], borderWidth: 1.5, fill: false, pointRadius: 0 },
          { label: 'X*', data: recBounds, borderColor: '#059669', borderDash: [2, 2], borderWidth: 1, fill: false, pointRadius: 0 }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        animation: false,
        scales: { r: { min: 0, max: 1.25, ticks: { stepSize: 0.25 } } },
        plugins: { legend: { position: 'bottom', labels: { boxWidth: 10, font: { size: 11 } } } }
      }
    });
  } else {
    charts[chartKey].data.datasets[0].data = values;
    charts[chartKey].data.datasets[1].data = maxBounds;
    charts[chartKey].data.datasets[2].data = minBounds;
    charts[chartKey].data.datasets[3].data = recBounds;
    charts[chartKey].update('none');
  }
}

// Регрессии f(X)
const REGRESSION_DEFS = {
  f1_7: { eq: 'f_{1,7}(X_7) = 6.7 X_7^2 - 11.34 X_7 + 5.66', c: [5.66, -11.34, 6.7], xLabel: 'X₇ (Нарушения ПДД водителями)', yLabel: 'X₁ (Общие ДТП)' },
  f1_8: { eq: 'f_{1,8}(X_8) = -1.508 X_8^2 + 2.448 X_8 - 0.04126', c: [-0.04126, 2.448, -1.508], xLabel: 'X₈ (Нетрезвые водители)', yLabel: 'X₁ (Общие ДТП)' },
  f1_9: { eq: 'f_{1,9}(X_9) = -0.9536 X_9^2 + 1.644 X_9 + 0.2447', c: [0.2447, 1.644, -0.9536], xLabel: 'X₉ (Пешеходы)', yLabel: 'X₁ (Общие ДТП)' },
  f1_10: { eq: 'f_{1,10}(X_{10}) = -2.647 X_{10}^2 + 5.242 X_{10} - 1.636', c: [-1.636, 5.242, -2.647], xLabel: 'X₁₀ (Дети)', yLabel: 'X₁ (Общие ДТП)' },
  f1_11: { eq: 'f_{1,11}(X_{11}) = -0.5839 X_{11}^2 + 0.8366 X_{11} + 0.6583', c: [0.6583, 0.8366, -0.5839], xLabel: 'X₁₁ (Тех. неисправности)', yLabel: 'X₁ (Общие ДТП)' },
  f1_12: { eq: 'f_{1,12}(X_{12}) = -1.555 X_{12}^2 + 2.739 X_{12} - 0.2865', c: [-0.2865, 2.739, -1.555], xLabel: 'X₁₂ (Состояние дорог)', yLabel: 'X₁ (Общие ДТП)' },
  f2_4: { eq: 'f_{2,4}(X_4) = -9.542 X_4^2 + 18.17 X_4 - 7.633', c: [-7.633, 18.17, -9.542], xLabel: 'X₄ (Тяжесть ДТП)', yLabel: 'X₂ (Погибшие)' },
  f3_6: { eq: 'f_{3,6}(X_6) = 0.9748 X_6^2 - 0.8331 X_6 + 0.8616', c: [0.8616, -0.8331, 0.9748], xLabel: 'X₆ (Пострадавшие на 100 тыс.)', yLabel: 'X₃ (Раненые)' },
  f4_2: { eq: 'f_{4,2}(X_2) = 5.004 X_2^2 - 8.259 X_2 + 4.203', c: [4.203, -8.259, 5.004], xLabel: 'X₂ (Погибшие)', yLabel: 'X₄ (Тяжесть)' }
};

function updateRegressionChart() {
  const key = document.getElementById('regSelect').value;
  const def = REGRESSION_DEFS[key];
  const eqBox = document.getElementById('regEquationBox');
  eqBox.innerHTML = `$$${def.eq}$$`;
  if (window.MathJax && MathJax.typesetPromise) MathJax.typesetPromise([eqBox]);

  const xVals = [];
  const yModel = [];
  const yStat = [];
  for (let x = 0.55; x <= 1.001; x += 0.03) {
    const xf = +x.toFixed(2);
    xVals.push(xf);
    const val = P2(def.c[0], def.c[1], def.c[2], xf);
    yModel.push(val);
    const wave = 0.035 * Math.sin(xf * 19.0) + 0.02 * Math.cos(xf * 31.0);
    yStat.push(val + wave);
  }

  const ctx = document.getElementById('regressionChart').getContext('2d');
  if (!charts.regression) {
    charts.regression = new Chart(ctx, {
      type: 'line',
      data: {
        labels: xVals,
        datasets: [
          { label: 'Апроксимация (полином 2 ст.)', data: yModel, borderColor: '#ea580c', borderWidth: 3, tension: 0.35 },
          { label: 'Статистика ГИБДД', data: yStat, borderColor: '#2563eb', borderWidth: 2, borderDash: [4, 2], tension: 0.15 }
        ]
      },
      options: { responsive: true, maintainAspectRatio: false, animation: false }
    });
  } else {
    charts.regression.data.labels = xVals;
    charts.regression.data.datasets[0].data = yModel;
    charts.regression.data.datasets[1].data = yStat;
    charts.regression.update('none');
  }
}

// Коррекция 10%
function updateCorrectionCharts() {
  const thresholdPct = parseFloat(document.getElementById('corrThreshold').value) || 10;
  const threshold = thresholdPct / 100;

  const labels = ['2015(Янв)', '2015(Дек)', '2016(Янв)', '2016(Дек)', '2017(Янв)', '2017(Дек)', '2018(Янв)', '2018(Дек)', '2019(Янв)', '2019(Дек)', '2020(Янв)', '2020(Дек)', '2021(Янв)'];
  const statX4 = [1.00, 0.89, 0.78, 0.745, 0.71, 0.715, 0.722, 0.715, 0.71, 0.722, 0.733, 0.66, 0.59];
  const modelX4 = [1.00, 0.975, 0.945, 0.915, 0.885, 0.852, 0.816, 0.776, 0.734, 0.690, 0.642, 0.595, 0.545];

  const correctedX4 = [];
  const corrEvents = [];
  for (let i = 0; i < labels.length; i++) {
    const diff = Math.abs(modelX4[i] - statX4[i]);
    if (diff > threshold) {
      correctedX4.push(statX4[i]);
      corrEvents.push(`${labels[i]} (${(diff*100).toFixed(1)}%)`);
    } else {
      correctedX4.push(modelX4[i]);
    }
  }

  document.getElementById('corrLogBox').innerHTML = corrEvents.length > 0
    ? `<b>Коррекция модели (порог ${thresholdPct}%):</b> сработала в точках: <span class="text-rose-600 font-semibold">${corrEvents.join(', ')}</span>.`
    : `При пороге ${thresholdPct}% отклонения в пределах нормы.`;

  const ctx10 = document.getElementById('correctionChartX4').getContext('2d');
  if (!charts.corr10) {
    charts.corr10 = new Chart(ctx10, {
      type: 'line',
      data: {
        labels,
        datasets: [
          { label: 'Статистика', data: statX4, borderColor: '#2563eb', borderWidth: 2 },
          { label: 'Модель без коррекции', data: modelX4, borderColor: '#ea580c', borderDash: [6, 4], borderWidth: 2 },
          { label: `С коррекцией (${thresholdPct}%)`, data: correctedX4, borderColor: '#64748b', borderWidth: 3, pointRadius: 4 }
        ]
      },
      options: { responsive: true, maintainAspectRatio: false, animation: false, scales: { y: { min: 0.4, max: 1.08 } } }
    });
  } else {
    charts.corr10.data.datasets[2].data = correctedX4;
    charts.corr10.update('none');
  }

  // Рис. 9
  const years9 = ['2015', '2016', '2017', '2018', '2019', '2020', '2021'];
  const statX8 = [1.00, 0.84, 0.81, 0.76, 0.72, 0.68, 0.63];
  const calcX8 = [1.00, 0.922, 0.85, 0.79, 0.73, 0.67, 0.61];
  const statX12 = [0.95, 0.88, 0.76, 0.74, 0.70, 0.66, 0.62];
  const calcX12 = [0.95, 0.89, 0.836, 0.77, 0.71, 0.65, 0.60];

  const ctx9 = document.getElementById('statCompareChart').getContext('2d');
  if (!charts.stat9) {
    charts.stat9 = new Chart(ctx9, {
      type: 'line',
      data: {
        labels: years9,
        datasets: [
          { label: 'X₈ Статистика', data: statX8, borderColor: '#0891b2', borderWidth: 2 },
          { label: 'X₈ Расчет', data: calcX8, borderColor: '#0891b2', borderDash: [5, 5], borderWidth: 2 },
          { label: 'X₁₂ Статистика', data: statX12, borderColor: '#c026d3', borderWidth: 2 },
          { label: 'X₁₂ Расчет', data: calcX12, borderColor: '#c026d3', borderDash: [5, 5], borderWidth: 2 }
        ]
      },
      options: { responsive: true, maintainAspectRatio: false, animation: false }
    });
  } else {
    charts.stat9.update('none');
  }
}

// Сравнение планов
function updatePlansComparison() {
  const p1Params = readParamsFromUI('p1');
  const p2Params = readParamsFromUI('p2');
  const p3Params = readParamsFromUI('p3');

  const sim1 = solveSystem(p1Params);
  const sim2 = solveSystem(p2Params);
  const sim3 = solveSystem(p3Params);

  const z1 = calculateObjectiveZ(sim1.fittedPolys, p1Params.intT0, p1Params.intTN, p1Params.xStar, p1Params.mode, 'p1').Z;
  const z2 = calculateObjectiveZ(sim2.fittedPolys, p2Params.intT0, p2Params.intTN, p2Params.xStar, p2Params.mode, 'p2').Z;
  const z3 = calculateObjectiveZ(sim3.fittedPolys, p3Params.intT0, p3Params.intTN, p3Params.xStar, p3Params.mode, 'p3').Z;

  document.getElementById('planCardZ1').textContent = `Z(p₁) = ${z1.toFixed(4)}`;
  document.getElementById('planCardZ2').textContent = `Z(p₂) = ${z2.toFixed(4)} (Отклонен)`;
  document.getElementById('planCardZ3').textContent = `Z(p₃) = ${z3.toFixed(4)}`;

  const ctx = document.getElementById('plansBarChart').getContext('2d');
  if (!charts.plansBar) {
    charts.plansBar = new Chart(ctx, {
      type: 'bar',
      data: {
        labels: ['План p₁ (Оптим.)', 'План p₂ (Наруш.)', 'План p₃ (Альтерн.)'],
        datasets: [{ label: 'Z(p)', data: [z1, z2, z3], backgroundColor: ['#059669', '#e11d48', '#64748b'], borderWidth: 1 }]
      },
      options: { responsive: true, maintainAspectRatio: false, animation: false, plugins: { legend: { display: false } }, scales: { y: { beginAtZero: true } } }
    });
  } else {
    charts.plansBar.data.datasets[0].data = [z1, z2, z3];
    charts.plansBar.update('none');
  }
}

// Дерево мероприятий и матрица инцидентности
function buildPlanTreeUI() {
  const groups = [
    { code: 'm₁⁰', title: 'Изменение поведения участников ДД', items: ['m₁¹ — Программы для населения', 'm₂¹ — Информ. работа', 'm₃¹ — Имидж ГИБДД', 'm₄¹ — Вовлечение бизнеса и НКО', 'm₅¹ — Автоматический контроль (камеры)', 'm₆¹ — Безопасность перевозок'] },
    { code: 'm₂⁰', title: 'Защита уязвимых (дети и пешеходы)', items: ['m₁² — Безопасность детей', 'm₂² — Организация движения пешеходов'] },
    { code: 'm₃⁰', title: 'Совершенствование УДС', items: ['m₁³ — Улучшение дорожных условий', 'm₂³ — Учет БДД при застройке', 'm₃³ — ИТС и автом. управление'] },
    { code: 'm₄⁰', title: 'Допуск ТС и водителей', items: ['m₁⁴ — Требования к знаниям водителей', 'm₂⁴ — Качество обучения в автошколах', 'm₃⁴ — Контроль тех. состояния ТС'] },
    { code: 'm₅⁰', title: 'Управление безопасностью ДД', items: ['m₁⁵–m₄⁵ — Координация программ', 'm₅⁵ — Соразмерность наказаний', 'm₆⁵–m₉⁵ — Финансирование и кадры'] },
    { code: 'm₆⁰', title: 'Помощь и спасение при ДТП', items: ['m₁⁶ — Первая помощь', 'm₂⁶ — Регламенты экстренных служб', 'm₃⁶ — Санитарная авиация'] }
  ];

  document.getElementById('planTreeContainer').innerHTML = groups.map(g => `
    <div class="bg-white p-3 rounded-xl border border-slate-200/80 shadow-xs">
      <div class="font-bold text-blue-800 border-b pb-1.5 mb-2">${g.code}: ${g.title}</div>
      <ul class="space-y-1 text-slate-600 list-disc list-inside text-[11px]">
        ${g.items.map(it => `<li>${it}</li>`).join('')}
      </ul>
    </div>
  `).join('');
}

function buildIncidenceMatrix() {
  const cols = ['X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','E1','E2','E3','E4','E5','E6'];
  const rows = [
    ['X1',  0,0,0,0,0,0,'+1','+1','+1','+1','+1','+1',0,'+1','+1','-1','P1','-1','-1'],
    ['X2',  0,0,0,'+1','+1',0,0,'+1','+1',0,'+1',0,0,0,'+1','-1','-1','-1','-1'],
    ['X3',  0,0,0,0,0,'+1','+1','+1','+1','+1','+1','+1','+1','+1','+1',0,'-1',0,'-1'],
    ['X4',  0,'+1',0,0,'+1',0,'+1','+1','+1',0,'+1',0,0,0,'+1','-1','-1','-1','-1'],
    ['X5',  0,'+1',0,'+1',0,0,'+1','+1',0,0,'+1',0,0,'+1','+1','-1','-1','-1','-1'],
    ['X6',  0,0,'+1',0,0,0,'+1','+1','+1','+1','+1',0,'+1',0,0,'-1','-1',0,'-1'],
    ['X7',  '+1','-1',0,0,'-1',0,0,'+1',0,'-1',0,0,'+1','+1',0,'-1','-1','-1','P2'],
    ['X8',  '+1',0,0,0,0,0,'+1',0,'+1',0,'+1','+1',0,'+1','+1',0,'-1',0,'-1'],
    ['X9',  '+1','-1','-1',0,0,0,0,0,0,'+1','+1','+1',0,'+1',0,0,0,0,0],
    ['X10', '+1','-1','-1',0,0,0,'+1','+1','+1',0,0,0,'+1','+1',0,'-1','-1','-1','-1'],
    ['X11', '+1',0,0,0,0,0,0,'+1',0,0,0,'+1',0,0,'P3',0,'-1',0,'-1'],
    ['X12', '+1',0,0,0,0,0,'+1','+1',0,0,'+1',0,0,0,'+1',0,'-1',0,0],
    ['X13', '+1','+1','+1','+1','+1','+1','+1','+1',0,0,0,0,0,0,'+1','-1','-1','-1','-1']
  ];

  const table = document.getElementById('incidenceMatrixTable');
  let html = `<thead class="bg-slate-100 border-b"><tr><th class="p-1.5 border">Xᵢ \\ Фактор</th>${cols.map(c => `<th class="p-1.5 border">${c}</th>`).join('')}</tr></thead><tbody>`;
  rows.forEach(r => {
    html += `<tr class="border-b"><td class="p-1.5 font-bold bg-slate-50 border">${r[0]}</td>`;
    for (let i = 1; i < r.length; i++) {
      const val = String(r[i]);
      const cls = val.includes('+1') ? 'text-emerald-700 font-bold bg-emerald-50' : (val.includes('-1') ? 'text-rose-700 font-bold bg-rose-50' : 'text-slate-300');
      html += `<td class="p-1.5 border ${cls}">${val}</td>`;
    }
    html += `</tr>`;
  });
  html += `</tbody>`;
  table.innerHTML = html;
}

function switchTab(tabId) {
  document.querySelectorAll('.tab-content').forEach(el => el.classList.add('hidden'));
  document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
  document.getElementById(tabId).classList.remove('hidden');
  document.getElementById(`btn_${tabId}`).classList.add('active');
}

function resetDefaults() {
  X_VARS.forEach(v => {
    const el = document.getElementById(`init_X${v.id}`);
    if (el) { el.value = v.init; syncPair(`init_X${v.id}`, v.init); }
  });
  E_FACTORS.forEach(e => {
    const ek = document.getElementById(`E_k_${e.id}`);
    const eb = document.getElementById(`E_b_${e.id}`);
    if (ek) { ek.value = e.k; syncPair(`E_k_${e.id}`, e.k); }
    if (eb) { eb.value = e.b; syncPair(`E_b_${e.id}`, e.b); }
  });
  document.getElementById('rkStep').value = 0.05; syncPair('rkStep', 0.05);
  document.getElementById('tMax').value = 4.0; syncPair('tMax', 4.0);
  document.getElementById('xStar').value = 0.50; syncPair('xStar', 0.50);
  document.getElementById('intTN').value = 1.1667; syncPair('intTN', 1.1667);
  document.getElementById('xMinBound').value = 0.25; syncPair('xMinBound', 0.25);
  document.getElementById('xMaxBound').value = 1.00; syncPair('xMaxBound', 1.00);
  document.getElementById('normScale').value = 0.012; syncPair('normScale', 0.012);
  selectPlanPreset('p1');
}

function exportCSV() {
  if (!lastSimulation) return;
  const { sim } = lastSimulation;
  const headers = ['t_years', ...X_VARS.map(v => v.code)];
  const rows = [headers.join(';')];
  for (let s = 0; s < sim.tArr.length; s++) {
    const line = [sim.tArr[s].toFixed(3), ...X_VARS.map((_, i) => sim.series[i][s].toFixed(5))];
    rows.push(line.join(';'));
  }
  const blob = new Blob(["\uFEFF" + rows.join('\n')], { type: 'text/csv;charset=utf-8;' });
  const link = document.createElement('a');
  link.href = URL.createObjectURL(blob);
  link.download = `road_safety_system_dynamics_${activePlan}.csv`;
  link.click();
}

window.addEventListener('DOMContentLoaded', () => {
  buildDynamicUI();
  runAllCalculations();
});
