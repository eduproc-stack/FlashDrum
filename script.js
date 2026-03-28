const els = {
  staticHeight: document.getElementById('staticHeight'),
  staticHeightNumber: document.getElementById('staticHeightNumber'),
  pipeDiameter: document.getElementById('pipeDiameter'),
  pipeDiameterNumber: document.getElementById('pipeDiameterNumber'),
  wheelDiameter: document.getElementById('wheelDiameter'),
  wheelDiameterNumber: document.getElementById('wheelDiameterNumber'),
  efficiency: document.getElementById('efficiency'),
  flowValue: document.getElementById('flowValue'),
  headValue: document.getElementById('headValue'),
  powerValue: document.getElementById('powerValue'),
  pumpEquation: document.getElementById('pumpEquation'),
  systemEquation: document.getElementById('systemEquation'),
  statusPill: document.getElementById('statusPill'),
  chart: document.getElementById('chart'),
  resetBtn: document.getElementById('resetBtn')
};

const defaults = { staticHeight: 25, pipeDiameter: 85, wheelDiameter: 200, efficiency: 60 };
const pumpBase = {
  a0_215: 47.03809524,
  a1_215: -0.0135873,
  a2_215: 0.00018667,
  a3_215: -0.00000676,
  a0_diff: 27.76349206,
  a1_diff: -0.00845503,
  a2_diff: 0.0001854,
  a3_diff: -0.00000913
};
const systemBase = {
  c0: [-0.55, -0.6119, -0.17762238, -0.06485],
  c1: [0.2667, 0.0725, 0.0191158, 0.007258],
  c2: [0.0227, 0.002825, 0.0006998, 0.000234],
  diameters: [50, 75, 100, 125],
  o0: { n34: -1.025, n35: 0.0076 },
  o1: { p0: 1.474705, p1: -0.0390597427, p2: 0.0003509472, p3: -0.0000010591 },
  o2: { p0: 0.180061, p1: -0.0053589827, p2: 0.0000528168, p3: -0.0000001716 }
};

function syncPair(rangeEl, numberEl) {
  const value = rangeEl.value;
  numberEl.value = value;
}
function syncReverse(rangeEl, numberEl) {
  const v = Number(numberEl.value);
  if (!Number.isFinite(v)) return;
  rangeEl.value = Math.min(Number(rangeEl.max), Math.max(Number(rangeEl.min), v));
  numberEl.value = rangeEl.value;
}

[['staticHeight','staticHeightNumber'],['pipeDiameter','pipeDiameterNumber'],['wheelDiameter','wheelDiameterNumber']].forEach(([a,b]) => {
  els[a].addEventListener('input', () => { syncPair(els[a], els[b]); update(); });
  els[b].addEventListener('input', () => { syncReverse(els[a], els[b]); update(); });
});
els.efficiency.addEventListener('input', update);
els.resetBtn.addEventListener('click', () => {
  Object.entries(defaults).forEach(([key, value]) => {
    if (els[key]) els[key].value = value;
    if (els[`${key}Number`]) els[`${key}Number`].value = value;
  });
  update();
});

function derivePumpCoefficients(wheelDiameter) {
  const ratio = (215 - wheelDiameter) / 50;
  return {
    a0: pumpBase.a0_215 - pumpBase.a0_diff * ratio,
    a1: pumpBase.a1_215 - pumpBase.a1_diff * ratio,
    a2: pumpBase.a2_215 - pumpBase.a2_diff * ratio,
    a3: pumpBase.a3_215 - pumpBase.a3_diff * ratio,
  };
}

function deriveSystemCoefficients(staticHeight, pipeDiameter) {
  const frictionA0 = Math.max(systemBase.o0.n34 + systemBase.o0.n35 * pipeDiameter, 0);
  const frictionA1 = systemBase.o1.p0 + systemBase.o1.p1 * pipeDiameter + systemBase.o1.p2 * pipeDiameter ** 2 + systemBase.o1.p3 * pipeDiameter ** 3;
  const frictionA2 = systemBase.o2.p0 + systemBase.o2.p1 * pipeDiameter + systemBase.o2.p2 * pipeDiameter ** 2 + systemBase.o2.p3 * pipeDiameter ** 3;
  return { b0: staticHeight + frictionA0, b1: frictionA1, b2: frictionA2 };
}

function pumpHead(q, c) { return c.a0 + c.a1 * q + c.a2 * q ** 2 + c.a3 * q ** 3; }
function systemHead(q, c) { return c.b0 + c.b1 * q + c.b2 * q ** 2; }

function solveOperationalPoint(pump, system) {
  const coeffs = [pump.a3, pump.a2 - system.b2, pump.a1 - system.b1, pump.a0 - system.b0];
  const roots = cubicRootsReal(coeffs);
  const valid = roots.filter((q) => q >= 0 && q <= 125);
  if (valid.length) {
    const q = valid.sort((a,b) => a-b)[0];
    return { flow: q, head: pumpHead(q, pump), inRange: true };
  }

  let best = { flow: 0, head: pumpHead(0, pump), diff: Infinity };
  for (let q = 0; q <= 125; q += 0.1) {
    const diff = Math.abs(pumpHead(q, pump) - systemHead(q, system));
    if (diff < best.diff) best = { flow: q, head: pumpHead(q, pump), diff };
  }
  return { flow: best.flow, head: best.head, inRange: best.diff < 1 };
}

function cubicRootsReal([a,b,c,d]) {
  if (Math.abs(a) < 1e-12) return quadraticRootsReal(b,c,d);
  const ba = b / a, ca = c / a, da = d / a;
  const p = ca - ba ** 2 / 3;
  const q = 2 * ba ** 3 / 27 - ba * ca / 3 + da;
  const disc = (q / 2) ** 2 + (p / 3) ** 3;
  if (disc > 1e-12) {
    const u = Math.cbrt(-q / 2 + Math.sqrt(disc));
    const v = Math.cbrt(-q / 2 - Math.sqrt(disc));
    return [u + v - ba / 3];
  }
  if (Math.abs(disc) <= 1e-12) {
    const u = Math.cbrt(-q / 2);
    return [2 * u - ba / 3, -u - ba / 3].filter((v, i, arr) => arr.findIndex((x) => Math.abs(x - v) < 1e-9) === i);
  }
  const r = Math.sqrt(-(p ** 3) / 27);
  const phi = Math.acos(-q / (2 * r));
  const m = 2 * Math.sqrt(-p / 3);
  return [0,1,2].map((k) => m * Math.cos((phi + 2 * Math.PI * k) / 3) - ba / 3);
}

function quadraticRootsReal(a,b,c) {
  if (Math.abs(a) < 1e-12) return Math.abs(b) < 1e-12 ? [] : [-c / b];
  const disc = b*b - 4*a*c;
  if (disc < 0) return [];
  if (disc === 0) return [-b / (2*a)];
  const s = Math.sqrt(disc);
  return [(-b + s)/(2*a), (-b - s)/(2*a)];
}

function fmt(value, digits = 2) {
  return Number(value).toLocaleString('sv-SE', { minimumFractionDigits: digits, maximumFractionDigits: digits });
}

function drawChart(pump, system, op) {
  const svg = els.chart;
  const width = 760, height = 440;
  const margin = { top: 24, right: 24, bottom: 52, left: 64 };
  const innerW = width - margin.left - margin.right;
  const innerH = height - margin.top - margin.bottom;
  const qs = Array.from({ length: 126 }, (_, i) => i);
  const data = qs.map((q) => ({ q, pump: pumpHead(q, pump), system: systemHead(q, system) }));
  const yMax = Math.max(5, ...data.flatMap((d) => [d.pump, d.system]), op.head) * 1.1;
  const x = (q) => margin.left + (q / 125) * innerW;
  const y = (h) => margin.top + innerH - (h / yMax) * innerH;
  const line = (key) => data.map((d, i) => `${i ? 'L' : 'M'}${x(d.q).toFixed(2)},${y(d[key]).toFixed(2)}`).join(' ');

  let grid = '';
  for (let i = 0; i <= 5; i++) {
    const gyVal = (yMax / 5) * i;
    const gy = y(gyVal);
    grid += `<line x1="${margin.left}" y1="${gy}" x2="${width - margin.right}" y2="${gy}" stroke="#d9e2ef" stroke-width="1" />`;
    grid += `<text x="${margin.left - 12}" y="${gy + 4}" text-anchor="end" font-size="12" fill="#5f7089">${fmt(gyVal,0)}</text>`;
  }
  for (let i = 0; i <= 5; i++) {
    const gxVal = 25 * i;
    const gx = x(gxVal);
    grid += `<line x1="${gx}" y1="${margin.top}" x2="${gx}" y2="${height - margin.bottom}" stroke="#eef3fa" stroke-width="1" />`;
    grid += `<text x="${gx}" y="${height - margin.bottom + 22}" text-anchor="middle" font-size="12" fill="#5f7089">${gxVal}</text>`;
  }

  svg.innerHTML = `
    <rect x="0" y="0" width="${width}" height="${height}" fill="transparent"></rect>
    ${grid}
    <line x1="${margin.left}" y1="${height - margin.bottom}" x2="${width - margin.right}" y2="${height - margin.bottom}" stroke="#10233f" stroke-width="1.5" />
    <line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${height - margin.bottom}" stroke="#10233f" stroke-width="1.5" />
    <text x="${width / 2}" y="${height - 14}" text-anchor="middle" font-size="13" fill="#5f7089">Flöde (m³/h)</text>
    <text x="18" y="${height / 2}" text-anchor="middle" font-size="13" fill="#5f7089" transform="rotate(-90 18 ${height / 2})">Lyfthöjd (mVP)</text>
    <path d="${line('pump')}" fill="none" stroke="#2357d8" stroke-width="3" stroke-linecap="round"></path>
    <path d="${line('system')}" fill="none" stroke="#14866d" stroke-width="3" stroke-linecap="round"></path>
    <circle cx="${x(op.flow)}" cy="${y(op.head)}" r="6" fill="#d13f5c"></circle>
    <text x="${Math.min(width - 90, x(op.flow) + 10)}" y="${Math.max(32, y(op.head) - 10)}" font-size="12" fill="#10233f">Driftpunkt ${fmt(op.flow,1)} / ${fmt(op.head,1)}</text>
  `;
}

function update() {
  const staticHeight = Number(els.staticHeight.value);
  const pipeDiameter = Number(els.pipeDiameter.value);
  const wheelDiameter = Number(els.wheelDiameter.value);
  const efficiency = Math.max(1, Number(els.efficiency.value) || defaults.efficiency) / 100;

  const pump = derivePumpCoefficients(wheelDiameter);
  const system = deriveSystemCoefficients(staticHeight, pipeDiameter);
  const op = solveOperationalPoint(pump, system);
  const power = (op.flow * op.head / 3600 * 9.81) / efficiency;
  const inArea = op.inRange && op.flow >= 0 && op.flow <= 125;

  els.flowValue.textContent = fmt(op.flow, 1);
  els.headValue.textContent = fmt(op.head, 1);
  els.powerValue.textContent = fmt(power, 2);
  els.pumpEquation.textContent = `H = ${fmt(pump.a0,4)} ${pump.a1 >= 0 ? '+' : '-'} ${fmt(Math.abs(pump.a1),4)}Q ${pump.a2 >= 0 ? '+' : '-'} ${fmt(Math.abs(pump.a2),6)}Q² ${pump.a3 >= 0 ? '+' : '-'} ${fmt(Math.abs(pump.a3),8)}Q³`;
  els.systemEquation.textContent = `H = ${fmt(system.b0,4)} ${system.b1 >= 0 ? '+' : '-'} ${fmt(Math.abs(system.b1),6)}Q ${system.b2 >= 0 ? '+' : '-'} ${fmt(Math.abs(system.b2),8)}Q²`;
  els.statusPill.textContent = inArea ? 'Inom område' : 'Kontrollera område';
  els.statusPill.className = `pill ${inArea ? 'ok' : 'warn'}`;
  drawChart(pump, system, op);
}

update();
