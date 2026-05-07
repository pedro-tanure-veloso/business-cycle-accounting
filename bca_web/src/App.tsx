import { useState, useEffect } from 'react';
import { 
  Activity, BarChart3, BrainCircuit, TrendingUp, TrendingDown, 
  AlertTriangle, Clock, ActivitySquare
} from 'lucide-react';
import { 
  BarChart, Bar, ComposedChart, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer,
  LineChart, Line, Legend
} from 'recharts';
import './App.css';

// Mock types for the JSON data contract
interface BCAData {
  quarter: string;
  macro_quarter: string;
  macro_overview: {
    gdp_growth_qoq: number;
    gdi_growth_qoq: number;
    supply_growth_qoq: number;
    gdp_growth_yoy: number;
    components: {
      consumption: { growth_qoq: number; contribution_to_gdp: number };
      investment: { growth_qoq: number; contribution_to_gdp: number };
      government: { growth_qoq: number; contribution_to_gdp: number };
      net_exports: { contribution_to_gdp: number };
    };
    historical_percentiles: {
      gdp: number;
      investment: number;
      consumption: number;
    };
  };
  time_series: {
    demand_contributions: Array<{ quarter: string; Consumption: number; Investment: number; Government: number; "Net Exports": number; "Total GDP Growth": number }>;
    supply_contributions: Array<{ quarter: string; Goods: number; Services: number; Government: number }>;
  };
  wedge_decomposition: {
    current_levels: Record<string, { sd_from_mean: number; percentile: number; trend: string }>;
    phi_statistics: Record<string, number>;
    historical_comparison: {
      most_analogous_episode: string;
      similarity_score: number;
    };
  };
  hypothesis_layer: {
    pattern_identification: string;
    candidate_mechanisms: Array<{
      wedge: string;
      mechanism: string;
      citations: string[];
      reasoning: string;
    }>;
    what_to_watch: Array<{
      mechanism: string;
      indicator: string;
    }>;
  };
}

function App() {
  const [data, setData] = useState<BCAData | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    // Fetch the production static data export
    fetch('/data/latest_quarter.json')
      .then(res => res.json())
      .then(jsonData => {
        setData(jsonData);
        setLoading(false);
      })
      .catch(err => {
        console.error("Failed to load data:", err);
        setLoading(false);
      });
  }, []);

  if (loading || !data) {
    return <div className="app-container" style={{ alignItems: 'center', justifyContent: 'center' }}>
      <ActivitySquare className="brand-icon" size={48} />
      <h2 style={{ marginTop: '1rem' }}>Loading BCA Data...</h2>
    </div>;
  }

  // Format phi-stats for the chart
  const phiChartData = Object.entries(data.wedge_decomposition.phi_statistics).map(([name, value]) => ({
    name: name.charAt(0).toUpperCase() + name.slice(1),
    value: value * 100 // convert to percentage
  }));

  const renderArrow = (value: number) => {
    return value >= 0 
      ? <TrendingUp size={28} color="var(--success)" style={{ marginLeft: '0.5rem' }} />
      : <TrendingDown size={28} color="var(--danger)" style={{ marginLeft: '0.5rem' }} />;
  };

  const annualizeRate = (val: number) => (Math.pow(1 + val, 4) - 1) * 100;
  const annualizeContrib = (val: number) => val * 4;

  return (
    <div className="app-container">
      {/* HEADER */}
      <header className="app-header">
        <div className="header-brand">
          <div className="brand-icon">
            <Activity size={24} />
          </div>
          <div>
            <h1 className="brand-title">Business Cycle Monitor</h1>
            <span className="brand-subtitle">US Quarterly Structural Snapshot</span>
          </div>
        </div>
        <div className="header-status">
          <div className="status-dot"></div>
          <span>Updated: {data.macro_quarter || data.quarter} (NIPA Release)</span>
        </div>
      </header>

      <main className="main-content">
        
        {/* SCREEN 1: MACRO OVERVIEW */}
        <section>
          <div className="section-header">
            <BarChart3 className="text-accent" size={24} color="var(--accent-primary)" />
            <h2 className="section-title">Macro Overview</h2>
            <div className="section-divider"></div>
          </div>
          
          <div className="grid-5">
            <div className="glass-panel kpi-card">
              <span className="kpi-label">GDP Growth</span>
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <span className="kpi-value">{annualizeRate(data.macro_overview.gdp_growth_qoq).toFixed(1)}%</span>
                {renderArrow(data.macro_overview.gdp_growth_qoq)}
              </div>
            </div>
            <div className="glass-panel kpi-card">
              <span className="kpi-label">Consumption</span>
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <span className="kpi-value">{annualizeRate(data.macro_overview.components.consumption.growth_qoq).toFixed(1)}%</span>
                {renderArrow(data.macro_overview.components.consumption.growth_qoq)}
              </div>
            </div>
            <div className="glass-panel kpi-card">
              <span className="kpi-label">Investment</span>
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <span className="kpi-value">{annualizeRate(data.macro_overview.components.investment.growth_qoq).toFixed(1)}%</span>
                {renderArrow(data.macro_overview.components.investment.growth_qoq)}
              </div>
            </div>
            <div className="glass-panel kpi-card">
              <span className="kpi-label">Government</span>
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <span className="kpi-value">{annualizeRate(data.macro_overview.components.government.growth_qoq).toFixed(1)}%</span>
                {renderArrow(data.macro_overview.components.government.growth_qoq)}
              </div>
            </div>
            <div className="glass-panel kpi-card">
              <span className="kpi-label">Net Exports</span>
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <span className="kpi-value">{annualizeContrib(data.macro_overview.components.net_exports.contribution_to_gdp).toFixed(1)}%</span>
                {renderArrow(data.macro_overview.components.net_exports.contribution_to_gdp)}
              </div>
            </div>
            <div className="glass-panel kpi-card">
              <span className="kpi-label">Income Optic (GDI)</span>
              <div style={{ display: 'flex', alignItems: 'center' }}>
                <span className="kpi-value">{annualizeRate(data.macro_overview.gdi_growth_qoq).toFixed(1)}%</span>
                {renderArrow(data.macro_overview.gdi_growth_qoq)}
              </div>
            </div>
          </div>
          <div style={{ marginTop: '1rem', fontSize: '0.875rem', color: 'var(--text-muted)' }}>
            * Note: All macro overview figures represent QoQ annualized rates or contributions to GDP.
          </div>

          <div style={{ marginTop: '2rem' }}>
            <div className="glass-panel">
              <h3 style={{ marginBottom: '1.5rem', color: 'var(--text-secondary)' }}>Demand Contributions to Growth</h3>
              <div className="chart-container">
                <ResponsiveContainer width="100%" height="100%">
                  <ComposedChart data={data.time_series.demand_contributions} stackOffset="sign" margin={{ top: 5, right: 5, left: -20, bottom: 5 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="var(--border-color)" vertical={false} />
                    <XAxis dataKey="quarter" stroke="var(--text-muted)" fontSize={12} />
                    <YAxis stroke="var(--text-muted)" fontSize={12} />
                    <Tooltip contentStyle={{ backgroundColor: 'var(--bg-surface)', borderColor: 'var(--border-color)' }} />
                    <Legend />
                    <Bar dataKey="Consumption" stackId="a" fill="#3b82f6" />
                    <Bar dataKey="Investment" stackId="a" fill="#8b5cf6" />
                    <Bar dataKey="Government" stackId="a" fill="#10b981" />
                    <Bar dataKey="Net Exports" stackId="a" fill="#f59e0b" />
                    <Line type="monotone" dataKey="Total GDP Growth" stroke="#ffffff" strokeWidth={2} dot={{ r: 3, fill: '#ffffff' }} />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>
            </div>
          </div>
        </section>

        {/* SCREEN 2: WEDGE DECOMPOSITION */}
        <section>
          <div className="section-header">
            <ActivitySquare className="text-accent" size={24} color="var(--accent-primary)" />
            <h2 className="section-title">Wedge Decomposition</h2>
            <div className="section-divider"></div>
          </div>
          <div style={{ color: 'var(--text-muted)', fontSize: '0.875rem', marginBottom: '1.5rem', fontStyle: 'italic' }}>
            * Note: Structural wedges are reported up to {data.quarter} due to auxiliary data release lags.
          </div>

          <div className="grid-2">
            <div className="glass-panel">
              <h3 style={{ marginBottom: '1.5rem', color: 'var(--text-secondary)' }}>f-Statistics: Output Explained by Wedge</h3>
              <div className="chart-container">
                <ResponsiveContainer width="100%" height="100%">
                  <BarChart data={phiChartData} layout="vertical" margin={{ top: 5, right: 30, left: 80, bottom: 5 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="var(--border-color)" horizontal={true} vertical={false} />
                    <XAxis type="number" unit="%" stroke="var(--text-muted)" />
                    <YAxis dataKey="name" type="category" stroke="var(--text-muted)" />
                    <Tooltip cursor={{ fill: 'rgba(255,255,255,0.05)' }} contentStyle={{ backgroundColor: 'var(--bg-surface)', borderColor: 'var(--border-color)', borderRadius: '8px' }} />
                    <Bar dataKey="value" fill="var(--accent-primary)" radius={[0, 4, 4, 0]} />
                  </BarChart>
                </ResponsiveContainer>
              </div>
            </div>

            <div className="glass-panel">
              <h3 style={{ marginBottom: '1.5rem', color: 'var(--text-secondary)' }}>Historical Context</h3>
              
              <div style={{ padding: '1rem', background: 'rgba(255,255,255,0.03)', borderRadius: 'var(--radius-md)', marginBottom: '1.5rem' }}>
                <span className="kpi-label" style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                  <Clock size={16} /> Most Analogous Historical Episode
                </span>
                <div style={{ display: 'flex', alignItems: 'flex-end', gap: '1rem', marginTop: '0.5rem' }}>
                  <span className="kpi-value" style={{ fontSize: '1.75rem' }}>
                    {data.wedge_decomposition.historical_comparison.most_analogous_episode}
                  </span>
                  <span className="badge badge-blue">
                    {(data.wedge_decomposition.historical_comparison.similarity_score * 100).toFixed(0)}% Match
                  </span>
                </div>
              </div>

              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem' }}>
                {Object.entries(data.wedge_decomposition.current_levels).map(([wedge, stats]) => (
                  <div key={wedge} style={{ padding: '1rem', border: '1px solid var(--border-color)', borderRadius: 'var(--radius-sm)' }}>
                    <div style={{ textTransform: 'capitalize', fontWeight: 500, marginBottom: '0.25rem' }}>{wedge} Wedge</div>
                    <div style={{ display: 'flex', justifyContent: 'space-between', color: 'var(--text-muted)', fontSize: '0.875rem' }}>
                      <span>{stats.sd_from_mean > 0 ? '+' : ''}{stats.sd_from_mean} SD</span>
                      <span>{stats.percentile}th pct</span>
                    </div>
                  </div>
                ))}
              </div>
            </div>
          </div>
        </section>

        {/* WEDGE EXPLANATION BOX */}
        <section>
          <div className="glass-panel" style={{ background: 'rgba(59, 130, 246, 0.05)', borderLeft: '4px solid var(--accent-primary)' }}>
            <h3 style={{ marginBottom: '1rem', display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
              <ActivitySquare size={20} color="var(--accent-primary)" /> What do these wedges represent?
            </h3>
            <div className="grid-2" style={{ gap: '1rem' }}>
              <div>
                <strong style={{ color: 'white' }}>Efficiency Wedge:</strong> Acts like a technology parameter. It captures anything that lowers aggregate productivity (e.g., input financing frictions, resource misallocation, or trade disruptions).
              </div>
              <div>
                <strong style={{ color: 'white' }}>Labor Wedge:</strong> Acts like a tax on labor income. It represents anything that drives a wedge between the marginal rate of substitution and the marginal product of labor (e.g., search frictions, sticky wages, labor taxes).
              </div>
              <div>
                <strong style={{ color: 'white' }}>Investment Wedge:</strong> Acts like a tax on investment. It captures anything that makes it harder or more expensive to turn today's resources into tomorrow's capital (distortions to the intertemporal Euler equation), such as firm financing frictions, liquidity constraints, or uncertainty shocks.
              </div>
              <div>
                <strong style={{ color: 'white' }}>Government Wedge:</strong> Represents output that is absorbed outside of domestic private consumption and investment. It is measured directly from the data as the sum of government consumption and net exports.
              </div>
            </div>
          </div>
        </section>

        {/* SCREEN 3: HYPOTHESIS LAYER */}
        <section>
          <div className="section-header">
            <BrainCircuit className="text-accent" size={24} color="var(--accent-secondary)" />
            <h2 className="section-title" style={{ color: 'var(--accent-secondary)' }}>Hypothesis Layer</h2>
            <div className="section-divider"></div>
          </div>

          <div className="glass-panel hypothesis-card">
            <div className="hypothesis-header">
              <Activity size={20} />
              <h3 style={{ margin: 0 }}>Pattern Identification</h3>
            </div>
            <p className="hypothesis-text">
              {data.hypothesis_layer.pattern_identification}
            </p>
          </div>

          <div className="grid-2" style={{ marginTop: '1.5rem' }}>
            <div className="glass-panel">
              <h3 style={{ marginBottom: '1.5rem', color: 'var(--text-secondary)' }}>Candidate Mechanisms</h3>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '1.5rem' }}>
                {data.hypothesis_layer.candidate_mechanisms.map((mech, idx) => (
                  <div key={idx} style={{ paddingBottom: '1.5rem', borderBottom: idx !== data.hypothesis_layer.candidate_mechanisms.length - 1 ? '1px solid var(--border-color)' : 'none' }}>
                    <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '0.5rem' }}>
                      <span className="badge badge-purple">{mech.wedge} Wedge</span>
                      <span style={{ fontSize: '0.75rem', color: 'var(--text-muted)' }}>{mech.citations.join(', ')}</span>
                    </div>
                    <h4 style={{ color: 'white' }}>{mech.mechanism}</h4>
                    <p style={{ color: 'var(--text-muted)', fontSize: '0.9rem', marginTop: '0.5rem' }}>{mech.reasoning}</p>
                  </div>
                ))}
              </div>
            </div>

            <div className="glass-panel">
              <h3 style={{ marginBottom: '1.5rem', color: 'var(--text-secondary)' }}>What to Watch</h3>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                {data.hypothesis_layer.what_to_watch.map((watch, idx) => (
                  <div key={idx} style={{ padding: '1rem', background: 'rgba(255,255,255,0.03)', borderRadius: 'var(--radius-md)' }}>
                    <h5 style={{ color: 'white', marginBottom: '0.25rem' }}>{watch.mechanism}</h5>
                    <p style={{ color: 'var(--text-muted)', fontSize: '0.9rem' }}>{watch.indicator}</p>
                  </div>
                ))}
              </div>
            </div>
          </div>

          <div className="ai-warning">
            <AlertTriangle size={20} />
            <strong>Hypotheses generated by AI — not structural findings.</strong> 
            Treat as a starting point for analysis, not a conclusion. Based on Gemini 3.1 Pro (High).
          </div>
        </section>

      </main>
    </div>
  );
}

export default App;
