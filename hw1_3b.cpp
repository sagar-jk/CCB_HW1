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
    return (g() >> 11) * (1.0 / 9007199254740992.0);
}

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
    vector<pair<int,int>> react;
    vector<pair<int,int>> delta;
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

    vector<string> tok;
    {
        istringstream iss(s);
        string t;
        while (iss >> t) tok.push_back(t);
    }
    for (size_t i = 0; i < tok.size(); i += 2)
        cnt[tok[i]] += stoi(tok[i+1]);

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

        vector<string> parts;
        string cur;
        for (char c : line) {
            if (c == ':') { parts.push_back(trim(cur)); cur.clear(); }
            else cur.push_back(c);
        }
        parts.push_back(trim(cur));

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
    for (int i = 0; i < (int)M.species.size(); i++)
        M.idx[M.species[i]] = i;

    for (auto& t : raw) {
        string lhs, rhs; double k;
        tie(lhs, rhs, k) = t;

        auto L = parse_pairs_side(lhs);
        auto R = parse_pairs_side(rhs);

        Reaction rx;
        rx.k = k;

        for (auto& kv : L)
            rx.react.push_back({M.idx[kv.first], kv.second});

        unordered_map<int,int> d;
        for (auto& kv : R) d[M.idx[kv.first]] += kv.second;
        for (auto& kv : L) d[M.idx[kv.first]] -= kv.second;

        for (auto& kv : d)
            if (kv.second != 0)
                rx.delta.push_back({kv.first, kv.second});

        M.rxns.push_back(rx);
    }

    return M;
}

static unordered_map<string,int> parse_lambda_in(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open " + path);

    unordered_map<string,int> init;
    string line;

    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        istringstream iss(line);
        string name; int val;
        iss >> name >> val;
        init[name] = val;
    }

    return init;
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

int main(int argc, char** argv) {

    if (argc < 3) {
        cout << "Usage: " << argv[0] << " problem.r problem.in\n";
        return 0;
    }

    Model M = parse_lambda_r(argv[1]);
    auto init = parse_lambda_in(argv[2]);

    vector<int> x(M.species.size(), 0);
    for (auto& kv : init)
        if (M.idx.count(kv.first))
            x[M.idx[kv.first]] = kv.second;

    if (!M.idx.count("Y")) {
        cerr << "Error: species 'Y' not found in reaction file.\n";
        return 1;
    }
    int Y_idx = M.idx["Y"];

    mt19937_64 gen(1);

    vector<double> a(M.rxns.size());

    double t = 0.0;
    double TMAX = 10000;

    while (t < TMAX) {

        double a0 = 0.0;
        for (int i = 0; i < (int)M.rxns.size(); i++) {
            a[i] = propensity(M.rxns[i], x);
            a0 += a[i];
        }

        if (a0 <= 0.0) break;

        double r = uni01(gen) * a0;
        double cum = 0.0;
        int chosen = 0;
        for (int i = 0; i < (int)a.size(); i++) {
            cum += a[i];
            if (cum >= r) { chosen = i; break; }
        }

        double u = uni01(gen);
        if (u <= 0.0) u = 1e-16;
        t += -log(u) / a0;

        for (auto [si, dv] : M.rxns[chosen].delta)
            x[si] += dv;
    }

    cout << "Final Y = " << x[Y_idx] << endl;
}