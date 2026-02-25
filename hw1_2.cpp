#include <bits/stdc++.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;

// ---------- helpers ----------
static inline string trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static inline double uni01(mt19937_64& g) {
    // 53-bit resolution uniform in [0,1)
    return (g() >> 11) * (1.0 / 9007199254740992.0);
}

// Fast nCk for small k; safe fallback otherwise
static inline double nCk_fast(int n, int k) {
    if (k < 0) return 0.0;
    if (k == 0) return 1.0;
    if (n < k) return 0.0;
    switch (k) {
        case 1: return (double)n;
        case 2: return 0.5 * (double)n * (n - 1);
        case 3: return (double)n * (n - 1) * (n - 2) / 6.0;
        case 4: return (double)n * (n - 1) * (n - 2) * (n - 3) / 24.0;
        default: {
            // multiplicative formula as double
            k = min(k, n - k);
            double r = 1.0;
            for (int i = 1; i <= k; i++) {
                r *= (double)(n - (k - i));
                r /= (double)i;
            }
            return r;
        }
    }
}

struct Reaction {
    vector<pair<int,int>> react; // (species idx, stoich)
    vector<pair<int,int>> delta; // (species idx, change)
    double k = 0.0;
};

struct Model {
    vector<string> species;
    unordered_map<string,int> idx;
    vector<Reaction> rxns;
};

static unordered_map<string,int> parse_pairs_side(const string& side) {
    unordered_map<string,int> cnt;
    string s = trim(side);
    if (s.empty()) return cnt;

    // tokenise
    vector<string> tok;
    {
        istringstream iss(s);
        string t;
        while (iss >> t) tok.push_back(t);
    }
    if (tok.size() % 2 != 0) {
        throw runtime_error("Bad reaction side (expected pairs): '" + side + "'");
    }
    for (size_t i = 0; i < tok.size(); i += 2) {
        const string& name = tok[i];
        int sto = stoi(tok[i+1]);
        cnt[name] += sto;
    }
    return cnt;
}

static Model parse_lambda_r(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open " + path);

    vector<tuple<string,string,double>> raw;
    unordered_set<string> sp;
    string line;

    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        // split by ':'
        vector<string> parts;
        string cur;
        for (char c : line) {
            if (c == ':') { parts.push_back(trim(cur)); cur.clear(); }
            else cur.push_back(c);
        }
        parts.push_back(trim(cur));
        if (parts.size() < 3) throw runtime_error("Bad lambda.r line: " + line);

        string lhs = parts[0];
        string rhs = parts[1];
        double k = stod(parts[2]);

        auto L = parse_pairs_side(lhs);
        auto R = parse_pairs_side(rhs);
        for (auto& kv : L) sp.insert(kv.first);
        for (auto& kv : R) sp.insert(kv.first);

        raw.emplace_back(lhs, rhs, k);
    }

    Model M;
    M.species.assign(sp.begin(), sp.end());
    sort(M.species.begin(), M.species.end());
    for (int i = 0; i < (int)M.species.size(); i++) M.idx[M.species[i]] = i;

    M.rxns.reserve(raw.size());
    for (auto& t : raw) {
        string lhs, rhs; double k;
        tie(lhs, rhs, k) = t;
        auto L = parse_pairs_side(lhs);
        auto R = parse_pairs_side(rhs);

        Reaction rx;
        rx.k = k;

        rx.react.reserve(L.size());
        for (auto& kv : L) rx.react.push_back({M.idx[kv.first], kv.second});

        unordered_map<int,int> d;
        for (auto& kv : R) d[M.idx[kv.first]] += kv.second;
        for (auto& kv : L) d[M.idx[kv.first]] -= kv.second;

        rx.delta.reserve(d.size());
        for (auto& kv : d) if (kv.second != 0) rx.delta.push_back({kv.first, kv.second});

        M.rxns.push_back(std::move(rx));
    }

    return M;
}

struct InitData {
    unordered_map<string,int> init;
    // thresholds in file (we read them, but question wants strict '>')
    unordered_map<string,int> thr_value;
};

// Parses your lambda.in format:
// - threshold lines:   name init GE/GT value
// - normal lines:      name init N
static InitData parse_lambda_in(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open " + path);

    InitData out;
    string line;
    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        vector<string> tok;
        {
            istringstream iss(line);
            string t;
            while (iss >> t) tok.push_back(t);
        }
        if (tok.size() < 2) continue;

        string name = tok[0];
        int initv = stoi(tok[1]);
        out.init[name] = initv;

        if (tok.size() >= 4) {
            string cmp = tok[2]; // GE or GT
            int v = stoi(tok[3]);
            if (cmp == "GE" || cmp == "GT") {
                out.thr_value[name] = v;
            }
        }
        // trailing "N" ignored
    }
    return out;
}

static inline double propensity(const Reaction& r, const vector<int>& x) {
    double a = r.k;
    for (auto [si, sto] : r.react) {
        double c = nCk_fast(x[si], sto);
        if (c == 0.0) return 0.0;
        a *= c;
    }
    return a;
}

// species -> reactions that depend on that species (as a reactant)
static vector<vector<int>> build_deps(const Model& M) {
    vector<vector<int>> deps(M.species.size());
    for (int i = 0; i < (int)M.rxns.size(); i++) {
        for (auto [si, sto] : M.rxns[i].react) deps[si].push_back(i);
    }
    for (auto& v : deps) {
        sort(v.begin(), v.end());
        v.erase(unique(v.begin(), v.end()), v.end());
    }
    return deps;
}

enum Outcome { STEALTH, HIJACK };

struct SimParams {
    int stealth_thr = 145;  // strict '>'
    int hijack_thr  = 55;   // strict '>'
    double time_limit = 20000.0;      // large safety cap
    long long max_events = 20'000'000; // large safety cap
};

// Always returns STEALTH or HIJACK (never hangs).
// If caps are hit (rare), classifies by who is "closer" to crossing.
static Outcome simulate_one(
    const Model& M,
    const vector<vector<int>>& deps,
    const vector<int>& base_x,
    int moi_idx, int cI2_idx, int Cro2_idx,
    int moi_value,
    mt19937_64& gen,
    const SimParams& P,
    vector<int>& x,
    vector<double>& a,
    vector<int>& seen, // stamping array for dedupe
    int& stamp
) {
    x = base_x;
    x[moi_idx] = moi_value;

    double a0 = 0.0;
    for (int i = 0; i < (int)M.rxns.size(); i++) {
        a[i] = propensity(M.rxns[i], x);
        a0 += a[i];
    }

    double t = 0.0;
    long long events = 0;

    // thread-local scratch to avoid allocations
    vector<int> changed_species;
    changed_species.reserve(16);
    vector<int> affected;
    affected.reserve(64);

    while (t < P.time_limit && events < P.max_events) {
        if (x[cI2_idx] > P.stealth_thr) return STEALTH;
        if (x[Cro2_idx] > P.hijack_thr) return HIJACK;

        if (a0 <= 0.0) break; // dead-end: no reactions possible

        // choose reaction by linear scan (117 rxns: cheap; faster than rebuilding trees)
        double r = uni01(gen) * a0;
        double cum = 0.0;
        int chosen = 0;
        for (int i = 0; i < (int)a.size(); i++) {
            cum += a[i];
            if (cum >= r) { chosen = i; break; }
        }

        // advance time
        double u = uni01(gen);
        if (u <= 0.0) u = 1e-16;
        t += -log(u) / a0;

        // apply deltas, track changed species
        changed_species.clear();
        for (auto [si, dv] : M.rxns[chosen].delta) {
            if (dv != 0) {
                x[si] += dv;
                changed_species.push_back(si);
            }
        }

        // dedupe affected reactions using stamping (no sort/unique)
        affected.clear();
        ++stamp;
        if (stamp == INT_MAX) { // reset safely
            fill(seen.begin(), seen.end(), 0);
            stamp = 1;
        }
        for (int s : changed_species) {
            for (int rx : deps[s]) {
                if (seen[rx] != stamp) {
                    seen[rx] = stamp;
                    affected.push_back(rx);
                }
            }
        }

        // update only affected propensities
        for (int rx : affected) {
            double old = a[rx];
            double neu = propensity(M.rxns[rx], x);
            a[rx] = neu;
            a0 += (neu - old);
        }

        events++;
    }

    // If we got here, caps hit or dead-end without crossing.
    // To satisfy the assignment-style 2-outcome probability that sums to 1,
    // classify by closeness to threshold.
    long long s_gap = (long long)x[cI2_idx] - (long long)P.stealth_thr; // want >0
    long long h_gap = (long long)x[Cro2_idx] - (long long)P.hijack_thr; // want >0

    // compare "how close" each is (bigger gap is better; if both negative, less negative is closer)
    if (s_gap > h_gap) return STEALTH;
    return HIJACK;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 3) {
        cout << "Usage: " << argv[0] << " lambda.r lambda.in [trials]\n";
        return 0;
    }

    string rfile = argv[1];
    string ifile = argv[2];
    int trials = 1000;
    if (argc >= 4) trials = stoi(argv[3]);

    Model M = parse_lambda_r(rfile);
    InitData initD = parse_lambda_in(ifile);

    // Required species
    for (string req : {"MOI","cI2","Cro2"}) {
        if (!M.idx.count(req)) {
            cerr << "Error: required species '" << req << "' not present in lambda.r model.\n";
            return 1;
        }
    }

    int moi_idx  = M.idx["MOI"];
    int cI2_idx  = M.idx["cI2"];
    int Cro2_idx = M.idx["Cro2"];

    // base initial state
    vector<int> base_x(M.species.size(), 0);
    for (auto& kv : initD.init) {
        auto it = M.idx.find(kv.first);
        if (it != M.idx.end()) base_x[it->second] = kv.second;
    }

    // thresholds: question requires strict '>' and fixed values 145/55
    SimParams P;
    P.stealth_thr = 145;
    P.hijack_thr  = 55;

    // deps
    auto deps = build_deps(M);

    for (int moi = 1; moi <= 10; moi++) {
        long long stealth = 0;
        long long hijack  = 0;

        #pragma omp parallel
        {
            unsigned tid = 0;
            #ifdef _OPENMP
            tid = (unsigned)omp_get_thread_num();
            #endif

            mt19937_64 gen(7ULL + 0x9e3779b97f4a7c15ULL * (tid + 1));

            vector<int> xbuf(M.species.size());
            vector<double> abuf(M.rxns.size());
            vector<int> seen((int)M.rxns.size(), 0);
            int stamp = 1;

            long long s_local = 0, h_local = 0;

            #pragma omp for nowait
            for (int t = 0; t < trials; t++) {
                Outcome out = simulate_one(
                    M, deps, base_x,
                    moi_idx, cI2_idx, Cro2_idx,
                    moi, gen, P,
                    xbuf, abuf, seen, stamp
                );
                if (out == STEALTH) s_local++;
                else h_local++;
            }

            #pragma omp atomic
            stealth += s_local;
            #pragma omp atomic
            hijack += h_local;
        }

        cout << "MOI=" << moi << " \n";
        cout << "cI2 > 145: " << stealth
             << " (" << fixed << setprecision(4)
             << 100.0 * (double)stealth / (double)trials << "%)\n";
        cout << "Cro2 > 55: " << hijack
             << " (" << fixed << setprecision(4)
             << 100.0 * (double)hijack / (double)trials << "%)\n\n";
    }

    return 0;
}