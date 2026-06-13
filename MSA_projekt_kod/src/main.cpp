#include "scoring.hpp"
#include "matrices.hpp"
#include "align.hpp"
#include "ga.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <filesystem>
#include <algorithm>

using namespace msa;
namespace fs = std::filesystem;

struct Instance { ScoreMode mode; std::vector<std::string> seqs; std::string name; };

static Instance readInstance(const std::string& path) {
    std::ifstream in(path);
    if (!in) { std::cerr << "Nie moge otworzyc " << path << "\n"; exit(1); }
    Instance inst; inst.name = fs::path(path).stem().string();
    std::string line; bool gotMode = false;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        while (!line.empty() && (line.back()=='\r'||line.back()=='\n'||line.back()==' ')) line.pop_back();
        if (line.empty()) continue;
        if (!gotMode) {
            inst.mode = (line == "PROT" || line == "PROTEIN") ? ScoreMode::PROTEIN : ScoreMode::DNA;
            gotMode = true;
        } else {
            inst.seqs.push_back(line);
        }
    }
    return inst;
}

static void writeInstance(const std::string& path, ScoreMode mode,
                          const std::vector<std::string>& seqs) {
    std::ofstream out(path);
    out << (mode == ScoreMode::PROTEIN ? "PROT" : "DNA") << "\n";
    for (auto& s : seqs) out << s << "\n";
}

static Scoring scoringFor(ScoreMode mode, int gopen, int gext) {
    if (mode == ScoreMode::PROTEIN) return makeBlosum62(gopen, gext);
    Scoring s; s.mode = ScoreMode::DNA;
    s.match = 1; s.mismatch = -1; s.gopen = gopen; s.gext = gext;
    return s;
}

static int cmdValidate() {
    int failures = 0;
    auto check = [&](const std::string& label, long got, long exp) {
        bool ok = (got == exp);
        std::cout << (ok ? "  [OK] " : "  [BLAD] ") << label
                  << " = " << got << " (oczekiwano " << exp << ")\n";
        if (!ok) ++failures;
    };

    std::cout << "== Walidacja 1: DNA, match=+1 mismatch=-1 gopen=-2 gext=-1 ==\n";
    {
        Scoring s; s.mode=ScoreMode::DNA; s.match=1; s.mismatch=-1; s.gopen=-2; s.gext=-1;
        std::vector<std::string> A = {
            "AGTCGTAG",
            "A-TCGTCG",
            "----GTAG",
            "-GTAG-AG"
        };
        check("w1-w2", scorePair(A[0],A[1],s),  3);
        check("w1-w3", scorePair(A[0],A[2],s), -1);
        check("w1-w4", scorePair(A[0],A[3],s),  0);
        check("w2-w3", scorePair(A[1],A[2],s), -3);
        check("w2-w4", scorePair(A[1],A[3],s), -5);
        check("w3-w4", scorePair(A[2],A[3],s), -3);
        check("SP calkowite", spScore(A,s), -9);
    }

    std::cout << "== Walidacja 2: bialka, uproszczony BLOSUM, prosta kara gap=-4 ==\n";
    {
        Scoring s = makeSimpleProtein(-4, -4);
        std::vector<std::string> A = {
            "AGTCGTAGNPST",
            "ASTCGTAG-PST",
            "-GTPG-AGN-ST",
            "A-NCGT--NP-T",
            "-G--ATA--PS-"
        };
        check("w1-w2", scorePair(A[0],A[1],s), 51);
        check("w1-w3", scorePair(A[0],A[2],s), 27);
        check("w1-w4", scorePair(A[0],A[3],s), 26);
        check("w1-w5", scorePair(A[0],A[4],s),  2);
        check("w3-w4", scorePair(A[2],A[3],s),-14);
        check("w4-w5", scorePair(A[3],A[4],s),-20);
        check("SP calkowite", spScore(A,s), 85);
    }

    std::cout << "== Walidacja 3: bialka, uproszczony BLOSUM, gopen=-12 gext=-1 ==\n";
    {
        Scoring s = makeSimpleProtein(-12, -1);
        std::vector<std::string> A = {
            "AGTCGTAGNPST",
            "ASTCGTAG-PST",
            "-GTPG-AGN-ST",
            "A-NCGT--NP-T",
            "-G--ATA--PS-"
        };
        check("w1-w2", scorePair(A[0],A[1],s), 43);
        check("w1-w3", scorePair(A[0],A[2],s),  3);
        check("w1-w4", scorePair(A[0],A[3],s),  5);
        check("w1-w5", scorePair(A[0],A[4],s),-24);
        check("SP calkowite", spScore(A,s), -216);
    }

    std::cout << (failures==0 ? "\nWSZYSTKIE TESTY ZALICZONE.\n"
                              : "\nUWAGA: sa bledy walidacji!\n");
    return failures == 0 ? 0 : 1;
}

// generator instancji
static std::string randSeq(int len, ScoreMode mode, RNG& rng) {
    const std::string dna = "ACGT";
    const std::string aa  = "ARNDCQEGHILKMFPSTWYV";
    const std::string& al = (mode==ScoreMode::PROTEIN) ? aa : dna;
    std::uniform_int_distribution<int> d(0, (int)al.size()-1);
    std::string s; s.reserve(len);
    for (int i = 0; i < len; ++i) s.push_back(al[d(rng)]);
    return s;
}

static std::string mutateSeq(const std::string& anc, ScoreMode mode,
                             double psub, double pind, int gaplen, RNG& rng) {
    const std::string dna = "ACGT";
    const std::string aa  = "ARNDCQEGHILKMFPSTWYV";
    const std::string& al = (mode==ScoreMode::PROTEIN) ? aa : dna;
    std::uniform_int_distribution<int> dch(0, (int)al.size()-1);
    std::uniform_real_distribution<double> u(0,1);
    std::string out;
    for (size_t i = 0; i < anc.size(); ++i) {
        if (u(rng) < pind) {
            int len = std::max(1, gaplen);
            if (u(rng) < 0.5) {
                i += (len - 1);
                continue;
            } else {
                for (int t = 0; t < len; ++t) out.push_back(al[dch(rng)]);
            }
        }
        char c = anc[i];
        if (u(rng) < psub) {
            char nc; do { nc = al[dch(rng)]; } while (nc == c && al.size() > 1);
            c = nc;
        }
        out.push_back(c);
    }
    if (out.empty()) out.push_back(al[dch(rng)]);
    return out;
}

// jedna instancja
static std::vector<std::string> genInstance(ScoreMode mode, int k, int ancLen,
                                            double psub, double pind, int gaplen,
                                            RNG& rng) {
    std::string anc = randSeq(ancLen, mode, rng);
    std::vector<std::string> seqs;
    for (int i = 0; i < k; ++i)
        seqs.push_back(mutateSeq(anc, mode, psub, pind, gaplen, rng));
    return seqs;
}

static int cmdGen(int argc, char** argv) {
    std::map<std::string,std::string> opt;
    for (int i = 2; i + 1 < argc; i += 2)
        if (std::string(argv[i]).rfind("--",0)==0) opt[std::string(argv[i]).substr(2)] = argv[i+1];

    auto geti = [&](const std::string&k,int d){ return opt.count(k)?std::stoi(opt[k]):d; };
    auto getd = [&](const std::string&k,double d){ return opt.count(k)?std::stod(opt[k]):d; };
    auto gets = [&](const std::string&k,std::string d){ return opt.count(k)?opt[k]:d; };

    std::string outdir = gets("outdir","instances");
    fs::create_directories(outdir);
    ScoreMode mode = (gets("mode","DNA")=="PROT") ? ScoreMode::PROTEIN : ScoreMode::DNA;
    unsigned seed = (unsigned)geti("seed", 12345);
    RNG rng(seed);

    int count = geti("count", 1);
    int k     = geti("k", 6);
    int len   = geti("len", 30);
    double psub = getd("psub", 0.10);
    double pind = getd("pind", 0.03);
    int gaplen  = geti("gaplen", 1);
    std::string prefix = gets("prefix","inst");

    for (int c = 0; c < count; ++c) {
        auto seqs = genInstance(mode, k, len, psub, pind, gaplen, rng);
        std::ostringstream nm;
        nm << outdir << "/" << prefix << "_" << c << ".txt";
        writeInstance(nm.str(), mode, seqs);
    }
    std::cout << "Wygenerowano " << count << " instancji w " << outdir
              << " (prefix=" << prefix << ", k=" << k << ", len=" << len
              << ", psub=" << psub << ", pind=" << pind << ", gaplen=" << gaplen << ")\n";
    return 0;
}

//uruchomienie na 1 instancji
static int cmdRun(int argc, char** argv) {
    if (argc < 3) { std::cerr << "run <plik> [--flagi]\n"; return 1; }
    Instance inst = readInstance(argv[2]);
    std::map<std::string,std::string> opt;
    for (int i = 3; i + 1 < argc; i += 2)
        if (std::string(argv[i]).rfind("--",0)==0) opt[std::string(argv[i]).substr(2)] = argv[i+1];
    auto geti=[&](const std::string&k,int d){return opt.count(k)?std::stoi(opt[k]):d;};
    auto getd=[&](const std::string&k,double d){return opt.count(k)?std::stod(opt[k]):d;};

    GAParams P;
    P.popSize=geti("pop",60); P.pCross=getd("pc",0.85); P.pMut=getd("pm",0.30);
    P.elitism=geti("elit",2); P.tournament=geti("tour",3); P.maxEvals=geti("evals",40000);
    int gopen=geti("gopen", inst.mode==ScoreMode::PROTEIN?-10:-2);
    int gext =geti("gext", -1);
    unsigned seed=(unsigned)geti("seed",1);

    Scoring sc = scoringFor(inst.mode, gopen, gext);
    RNG rng(seed);
    GAResult R = runGA(inst.seqs, sc, P, rng);

    std::cout << "Instancja: " << inst.name << "  k=" << inst.seqs.size()
              << "  tryb=" << (inst.mode==ScoreMode::PROTEIN?"PROT":"DNA") << "\n";
    std::cout << "Init best SP = " << R.initBest
              << "   Final best SP = " << R.finalBest
              << "   (poprawa " << (R.finalBest - R.initBest) << ")\n";
    std::cout << "Generacje=" << R.gens << "  oceny=" << R.evals << "\n";
    std::cout << "Najlepsze dopasowanie (dlugosc " << alnLen(R.bestAlign) << "):\n";
    for (auto& r : R.bestAlign) std::cout << "  " << r << "\n";
    std::cout << "Poprawne? " << (isValid(R.bestAlign, inst.seqs) ? "TAK" : "NIE") << "\n";
    return 0;
}

static std::vector<std::string> collect(const std::string& dir, const std::string& prefix) {
    std::vector<std::string> out;
    for (auto& e : fs::directory_iterator(dir)) {
        std::string fn = e.path().filename().string();
        if (fn.rfind(prefix,0)==0 && e.path().extension()==".txt") out.push_back(e.path().string());
    }
    std::sort(out.begin(), out.end());
    return out;
}

//parametry
static int cmdSweep(int argc, char** argv) {
    if (argc < 4) { std::cerr << "sweep <param> <out.csv> [--flagi]\n"; return 1; }
    std::string param = argv[2];
    std::string outcsv = argv[3];
    std::map<std::string,std::string> opt;
    for (int i = 4; i + 1 < argc; i += 2)
        if (std::string(argv[i]).rfind("--",0)==0) opt[std::string(argv[i]).substr(2)] = argv[i+1];
    auto geti=[&](const std::string&k,int d){return opt.count(k)?std::stoi(opt[k]):d;};
    auto getd=[&](const std::string&k,double d){return opt.count(k)?std::stod(opt[k]):d;};
    auto gets=[&](const std::string&k,std::string d){return opt.count(k)?opt[k]:d;};

    std::string dir = gets("dir","instances");
    std::string prefix = gets("prefix","small");
    int runs = geti("runs",10);
    int gopen = geti("gopen", -2), gext = geti("gext", -1);

    GAParams base;
    base.popSize=geti("pop",60); base.pCross=getd("pc",0.85); base.pMut=getd("pm",0.30);
    base.elitism=geti("elit",2); base.tournament=geti("tour",3); base.maxEvals=geti("evals",40000);

    std::vector<std::string> vals;
    { std::stringstream ss(gets("values","")); std::string tok;
      while (std::getline(ss,tok,',')) if(!tok.empty()) vals.push_back(tok); }

    auto files = collect(dir, prefix);
    if (files.empty()) { std::cerr << "Brak instancji "<<prefix<<" w "<<dir<<"\n"; return 1; }

    std::ofstream csv(outcsv);
    csv << "param,value,instance,run,init_best,final_best,improvement,gens,evals\n";

    for (auto& vs : vals) {
        for (auto& f : files) {
            Instance inst = readInstance(f);
            Scoring sc = scoringFor(inst.mode, gopen, gext);
            for (int r = 0; r < runs; ++r) {
                GAParams P = base;
                if      (param=="pop")  P.popSize = std::stoi(vs);
                else if (param=="pm")   P.pMut    = std::stod(vs);
                else if (param=="pc")   P.pCross  = std::stod(vs);
                else if (param=="elit") P.elitism = std::stoi(vs);
                else if (param=="tour") P.tournament = std::stoi(vs);
                else { std::cerr<<"Nieznany param "<<param<<"\n"; return 1; }
                RNG rng((unsigned)(1000*r + 7));
                GAResult R = runGA(inst.seqs, sc, P, rng);
                csv << param << "," << vs << "," << inst.name << "," << r << ","
                    << R.initBest << "," << R.finalBest << ","
                    << (R.finalBest-R.initBest) << "," << R.gens << "," << R.evals << "\n";
            }
        }
        std::cerr << "  ..."<<param<<"="<<vs<<" gotowe\n";
    }
    std::cout << "Zapisano " << outcsv << "\n";
    return 0;
}

//różne klasy
static int cmdClasses(int argc, char** argv) {
    if (argc < 3) { std::cerr << "classes <out.csv> [--flagi]\n"; return 1; }
    std::string outcsv = argv[2];
    std::map<std::string,std::string> opt;
    for (int i = 3; i + 1 < argc; i += 2)
        if (std::string(argv[i]).rfind("--",0)==0) opt[std::string(argv[i]).substr(2)] = argv[i+1];
    auto geti=[&](const std::string&k,int d){return opt.count(k)?std::stoi(opt[k]):d;};
    auto getd=[&](const std::string&k,double d){return opt.count(k)?std::stod(opt[k]):d;};
    auto gets=[&](const std::string&k,std::string d){return opt.count(k)?opt[k]:d;};

    std::string dir = gets("dir","instances");
    int runs = geti("runs",10);
    int gopen = geti("gopen", -2), gext = geti("gext", -1);
    GAParams P; P.popSize=geti("pop",60); P.pCross=getd("pc",0.85);
    P.pMut=getd("pm",0.30); P.elitism=geti("elit",2); P.tournament=geti("tour",3);
    P.maxEvals=geti("evals",40000);

    std::vector<std::string> classes;
    { std::stringstream ss(gets("classes","")); std::string t;
      while(std::getline(ss,t,',')) if(!t.empty()) classes.push_back(t); }

    std::ofstream csv(outcsv);
    csv << "class,instance,run,init_best,final_best,improvement,evals\n";
    for (auto& cl : classes) {
        auto files = collect(dir, cl);
        for (auto& f : files) {
            Instance inst = readInstance(f);
            Scoring sc = scoringFor(inst.mode, gopen, gext);
            for (int r=0;r<runs;++r){
                RNG rng((unsigned)(1000*r+7));
                GAResult R = runGA(inst.seqs, sc, P, rng);
                csv << cl << "," << inst.name << "," << r << ","
                    << R.initBest << "," << R.finalBest << ","
                    << (R.finalBest-R.initBest) << "," << R.evals << "\n";
            }
        }
        std::cerr << "  ...klasa "<<cl<<" gotowa\n";
    }
    std::cout << "Zapisano " << outcsv << "\n";
    return 0;
}

//zbieżność
static int cmdConverge(int argc, char** argv) {
    if (argc < 3) { std::cerr << "converge <out.csv> [--flagi]\n"; return 1; }
    std::string outcsv = argv[2];
    std::map<std::string,std::string> opt;
    for (int i = 3; i + 1 < argc; i += 2)
        if (std::string(argv[i]).rfind("--",0)==0) opt[std::string(argv[i]).substr(2)] = argv[i+1];
    auto geti=[&](const std::string&k,int d){return opt.count(k)?std::stoi(opt[k]):d;};
    auto getd=[&](const std::string&k,double d){return opt.count(k)?std::stod(opt[k]):d;};
    auto gets=[&](const std::string&k,std::string d){return opt.count(k)?opt[k]:d;};

    std::string dir = gets("dir","instances");
    std::string prefix = gets("prefix","small");
    int runs = geti("runs",10);
    int gopen = geti("gopen",-2), gext=geti("gext",-1);
    GAParams P; P.popSize=geti("pop",60); P.pCross=getd("pc",0.85);
    P.pMut=getd("pm",0.30); P.elitism=geti("elit",2); P.tournament=geti("tour",3);
    P.maxEvals=geti("evals",40000);

    auto files = collect(dir, prefix);
    int pts = geti("points",40);
    std::vector<long> grid(pts);
    for (int i=0;i<pts;++i) grid[i] = (long)((double)(i+1)/pts * P.maxEvals);

    std::vector<std::vector<double>> samples(pts);
    for (auto& f : files) {
        Instance inst = readInstance(f);
        Scoring sc = scoringFor(inst.mode, gopen, gext);
        for (int r=0;r<runs;++r){
            RNG rng((unsigned)(1000*r+7));
            GAResult R = runGA(inst.seqs, sc, P, rng);
            int gi = 0;
            for (int i=0;i<pts;++i){
                while (gi+1 < (int)R.histEvals.size() && R.histEvals[gi+1] <= grid[i]) ++gi;
                double v = R.histBest.empty()? (double)R.initBest : (double)R.histBest[std::min((int)R.histBest.size()-1,gi)];
                samples[i].push_back(v);
            }
        }
    }
    std::ofstream csv(outcsv);
    csv << "evals,best_mean,best_std\n";
    for (int i=0;i<pts;++i){
        double m=0; for(double v:samples[i]) m+=v; m/=samples[i].size();
        double s=0; for(double v:samples[i]) s+=(v-m)*(v-m); s=std::sqrt(s/samples[i].size());
        csv << grid[i] << "," << m << "," << s << "\n";
    }
    std::cout << "Zapisano " << outcsv << "\n";
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Uzycie: " << argv[0]
                  << " {validate|gen|run|sweep|classes|converge} ...\n";
        return 1;
    }
    std::string cmd = argv[1];
    if (cmd=="validate") return cmdValidate();
    if (cmd=="gen")      return cmdGen(argc, argv);
    if (cmd=="run")      return cmdRun(argc, argv);
    if (cmd=="sweep")    return cmdSweep(argc, argv);
    if (cmd=="classes")  return cmdClasses(argc, argv);
    if (cmd=="converge") return cmdConverge(argc, argv);
    std::cerr << "Nieznane polecenie: " << cmd << "\n";
    return 1;
}
