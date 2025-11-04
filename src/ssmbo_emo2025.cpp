
#include <iostream>
#include <filesystem>
#include <tigon/Tigon.h>
extern "C"
{
#include <coco/coco.h>
}

using namespace Tigon;
using namespace Tigon::Representation;
using namespace Tigon::Operators;
using namespace Tigon::Algorithms;
using namespace Tigon::Log;

#define LOG_POPULATION

inline void dispVector(TVector<double> vec, string sep = "\t", string endLine = "\n")
{
    for(size_t i=0; i<vec.size(); i++) {
        cout << vec[i];
        if(i<vec.size()-1) {
            cout << sep;
        }
    }
    cout << endLine;
}

IFunctionSPtr BiObjProblem(int functionID=1, int nInputs=2, int instanceID=1,
                           bool useObserver=false, TString algorithmName="")
{
    BBOBBiObj* function = new BBOBBiObj();
    function->TP_defineFunctionID(functionID);
    function->TP_defineInstanceID(instanceID);
    function->TP_defineNInputs(nInputs);
    function->TP_defineUseObserver(useObserver);
    function->TP_defineObserverAlgorithmName(algorithmName);

    IFunctionSPtr func(function);
    return func;
}

/**
 * @brief Determines the smallest (normalized) Euclidean distance between all
 * solutions in solution set and the region of interest (RI)
 */
double minDistRI(const TVector<IMappingSPtr>& solutions,
                 const TVector<double>& ideal,
                 const TVector<double>& nadir)
{
    double minDist = 0.0;
    bool firstValue = true;
    for(const IMappingSPtr& sol : solutions) {
        TVector<double> f = sol->doubleObjectiveVec();

        TVector<double> ri(f.size());  // region of interest
        for(size_t i=0; i<f.size(); i++) {
            if(f[i] - ri[i] > 0) {
                ri[i] = nadir[i];
            } else {
                ri[i] = f[i];
            }
        }

        double dist = 0.0;
        for(size_t i=0; i<f.size(); i++) {
            dist += std::pow(fabs((f[i] - ri[i]) / (nadir[i] - ideal[i])), 2.0);
        }
        dist = std::sqrt(dist);

        if(firstValue) {
            minDist = dist;
            firstValue = false;
        } else {
            minDist = std::min(dist, minDist);
        }
    }
    return minDist;
}

double hypervolumeTigon(const ISet* archive,
                        const TVector<double>& ideal,
                        const TVector<double>& nadir)
{
    TVector<TVector<double>> TArchive;
    for(const IMappingSPtr& sol : archive->all()) {
        TVector<double> z = sol->doubleObjectiveVec();
        TVector<double> znorm(z.size());
        for(size_t i=0; i<z.size(); i++) {
            znorm[i] = (z[i] - ideal[i])/(nadir[i] - ideal[i]);
        }
        TArchive.push_back(znorm);
    }
    return hypervolume(TArchive, TVector<double>(2, 1.0));
}

void bmk_Random_BBOBBIOBJ(int budget, int seed, int functionID, int nInputs,
                         int instanceID)
{
    TRAND.defineSeed(seed);
    cout.precision(5);

    PSetBase*     base      = new PSetBase();
    IFormulation* form      = new IFormulation(base);
    RandomInit*    init     = new RandomInit(form);
    Evaluator*    evaluator = new Evaluator(init);
    ITermination* tend      = new ITermination(evaluator);

    base->defineKeepArchive(true);

    bool useObserver = false;
    IFunctionSPtr func = BiObjProblem(functionID, nInputs, instanceID,
                                      useObserver, "Random");
    ProblemSPtr prob(new Problem);
    prob->appendFunction(func);
    prob->processProblemDefinition();

    TVector<double> nadir = static_cast<BBOBBiObj*>(func.get())->nadirVector();
    TVector<double> ideal = static_cast<BBOBBiObj*>(func.get())->idealVector();
    double bestHyp = static_cast<BBOBBiObj*>(func.get())->bestHypervolumeValue();

    form->defineProblem(prob);

    init->TP_defineSetSize(budget);

    tend->defineBudget(1);

    LogManagerSPtr log = tend->log();
    log->defineLogAll(true);

    coco_archive_t* arc = coco_archive("bbob-biobj", functionID, nInputs, instanceID);

    TVector<double> Icoco;
    TVector<double> IcocoComp;

    tend->evaluate();
    ISet* evaluatedMappings = init->outputSet(0);
    TVector<IMappingSPtr> currentEvalMappings;
    for(int i=0; i<evaluatedMappings->size(); i++) {
        IMappingSPtr sol = evaluatedMappings->at(i);
        TString sName = "s" + std::to_string(i);
        TVector<double> z = sol->doubleObjectiveVec();
        coco_archive_add_solution(arc, z[0], z[1], sName.data());
        currentEvalMappings.push_back(sol);

        int nSolsArchive = coco_archive_get_number_of_solutions(arc);
        if(nSolsArchive == 2) {
            double minDist = minDistRI(currentEvalMappings, ideal, nadir);
            Icoco.push_back(minDist);
            IcocoComp.push_back(minDist);
        } else {
            double hpv = -coco_archive_get_hypervolume(arc);
            Icoco.push_back(hpv);
            ISet* archive = base->setWithTag(Tigon::TNonDominatedArchive);
            hpv = -hypervolumeTigon(archive, ideal, nadir);
            IcocoComp.push_back(hpv);
        }
    }

#ifdef LOG_POPULATION
    // Save population and evaluations
    JsonObject customFields;
    customFields["Random Seed"] = seed;

    JsonArray nadirArray = {nadir[0], nadir[1]};
    customFields["Nadir Vector"] = nadirArray;

    JsonArray idealArray = {ideal[0], ideal[1]};
    customFields["Ideal Vector"] = idealArray;

    customFields["Best Hypervolume"] = bestHyp;

    JsonArray hypervolumeArray;
    for(double hpv : Icoco) {
        hypervolumeArray.push_back(hpv);
    }
    customFields["ICoco"] = hypervolumeArray;

    JsonArray hypervolumeArrayComp;
    for(double hpv : IcocoComp) {
        hypervolumeArrayComp.push_back(hpv);
    }
    customFields["ICocoComp"] = hypervolumeArrayComp;

    log->logPopulation(evaluatedMappings, "Main Population", customFields);

    TString filePath = "outputs/Random/D" + std::to_string(nInputs) + "/runD" + std::to_string(nInputs) + "F" + std::to_string(functionID);
    std::filesystem::create_directories(filePath);

    JsonDocument jdoc;
    jdoc.setObject(log->populationsLog());
    TString fileName;
    fileName = filePath + "/population_Random_BBOB-BIOBJ_F"
               + std::to_string(functionID) + "_D"
               + std::to_string(nInputs) + "_I"
               + std::to_string(instanceID) + "_S"
               + std::to_string(seed) + ".json";
    Tigon::writeToFile(TString(jdoc.toJson().c_str()), fileName);

    JsonDocument jdoc2;
    jdoc2.setObject(log->evaluationsLog());
    TString fileName2;
    fileName2 = filePath + "/evaluations_Random_BBOB-BIOBJ_F"
                + std::to_string(functionID) + "_D"
                + std::to_string(nInputs) + "_I"
                + std::to_string(instanceID) + "_S"
                + std::to_string(seed) + ".json";
    Tigon::writeToFile(TString(jdoc2.toJson().c_str()), fileName2);
#endif

    coco_archive_free(arc);

    delete tend;
    delete evaluator;
    delete init;
    delete form;
    delete base;
}

void bmk_Sobol_BBOBBIOBJ(int budget, int seed, int functionID, int nInputs,
                         int instanceID)
{
    TRAND.defineSeed(seed);
    cout.precision(5);

    PSetBase*     base      = new PSetBase();
    IFormulation* form      = new IFormulation(base);
    SobolInit*    init      = new SobolInit(form);
    Evaluator*    evaluator = new Evaluator(init);
    ITermination* tend      = new ITermination(evaluator);

    base->defineKeepArchive(true);

    bool useObserver = false;
    IFunctionSPtr func = BiObjProblem(functionID, nInputs, instanceID,
                                      useObserver, "Sobol");
    ProblemSPtr prob(new Problem);
    prob->appendFunction(func);
    prob->processProblemDefinition();

    TVector<double> nadir = static_cast<BBOBBiObj*>(func.get())->nadirVector();
    TVector<double> ideal = static_cast<BBOBBiObj*>(func.get())->idealVector();
    double bestHyp = static_cast<BBOBBiObj*>(func.get())->bestHypervolumeValue();

    form->defineProblem(prob);

    init->TP_defineSetSize(budget);

    tend->defineBudget(1);

    LogManagerSPtr log = tend->log();
    log->defineLogAll(true);

    coco_archive_t* arc = coco_archive("bbob-biobj", functionID, nInputs, instanceID);

    TVector<double> Icoco;
    TVector<double> IcocoComp;

    tend->evaluate();
    ISet* evaluatedMappings = init->outputSet(0);
    TVector<IMappingSPtr> currentEvalMappings;
    for(int i=0; i<evaluatedMappings->size(); i++) {
        IMappingSPtr sol = evaluatedMappings->at(i);
        TString sName = "s" + std::to_string(i);
        TVector<double> z = sol->doubleObjectiveVec();
        coco_archive_add_solution(arc, z[0], z[1], sName.data());
        currentEvalMappings.push_back(sol);

        int nSolsArchive = coco_archive_get_number_of_solutions(arc);
        if(nSolsArchive == 2) {
            double minDist = minDistRI(currentEvalMappings, ideal, nadir);
            Icoco.push_back(minDist);
            IcocoComp.push_back(minDist);
        } else {
            double hpv = -coco_archive_get_hypervolume(arc);
            Icoco.push_back(hpv);
            ISet* archive = base->setWithTag(Tigon::TNonDominatedArchive);
            hpv = -hypervolumeTigon(archive, ideal, nadir);
            IcocoComp.push_back(hpv);
        }
    }

#ifdef LOG_POPULATION
    // Save population and evaluations
    JsonObject customFields;
    customFields["Random Seed"] = seed;

    JsonArray nadirArray = {nadir[0], nadir[1]};
    customFields["Nadir Vector"] = nadirArray;

    JsonArray idealArray = {ideal[0], ideal[1]};
    customFields["Ideal Vector"] = idealArray;

    customFields["Best Hypervolume"] = bestHyp;

    JsonArray hypervolumeArray;
    for(double hpv : Icoco) {
        hypervolumeArray.push_back(hpv);
    }
    customFields["ICoco"] = hypervolumeArray;

    JsonArray hypervolumeArrayComp;
    for(double hpv : IcocoComp) {
        hypervolumeArrayComp.push_back(hpv);
    }
    customFields["ICocoComp"] = hypervolumeArrayComp;

    log->logPopulation(evaluatedMappings, "Main Population", customFields);

    TString filePath = "outputs/Sobol/D" + std::to_string(nInputs) + "/runD" + std::to_string(nInputs) + "F" + std::to_string(functionID);
    std::filesystem::create_directories(filePath);

    JsonDocument jdoc;
    jdoc.setObject(log->populationsLog());
    TString fileName;
    fileName = filePath + "outputs/population_Sobol_BBOB-BIOBJ_F"
               + std::to_string(functionID) + "_D"
               + std::to_string(nInputs) + "_I"
               + std::to_string(instanceID) + "_S"
               + std::to_string(seed) + ".json";
    Tigon::writeToFile(TString(jdoc.toJson().c_str()), fileName);

    JsonDocument jdoc2;
    jdoc2.setObject(log->evaluationsLog());
    TString fileName2;
    fileName2 = filePath + "outputs/evaluations_Sobol_BBOB-BIOBJ_F"
                + std::to_string(functionID) + "_D"
                + std::to_string(nInputs) + "_I"
                + std::to_string(instanceID) + "_S"
                + std::to_string(seed) + ".json";
    Tigon::writeToFile(TString(jdoc2.toJson().c_str()), fileName2);
#endif

    coco_archive_free(arc);

    delete tend;
    delete evaluator;
    delete init;
    delete form;
    delete base;
}

void bmk_ParEGO_BBOBBIOBJ(int budget, int seed, int functionID, int nInputs, int instanceID,
                          Tigon::ScalarisationType sFunc=Tigon::WeightedChebyshevAugmented)
{
    TRAND.defineSeed(seed);
    cout.precision(5);

    int   nDirs    = 11;
    int   popSize  = 11 * nInputs - 1;
    // int   maxSurrogateSize = popSize;
    int   maxSurrogateSize = 100;

    PSetBase*     base      = new PSetBase();
    IFormulation* form      = new IFormulation(base);
    SobolInit*    init      = new SobolInit(form);
    Evaluator*    evaluator = new Evaluator(init);
    ParEGO*       alg       = new ParEGO(evaluator);
    ITermination* tend      = new ITermination(alg);

    base->defineKeepArchive(true);

    bool useObserver = false;
    IFunctionSPtr func = BiObjProblem(functionID, nInputs, instanceID,
                                      useObserver, "ParEGO");
    ProblemSPtr prob(new Problem);
    prob->appendFunction(func);
    prob->processProblemDefinition();

    TVector<double> nadir = static_cast<BBOBBiObj*>(func.get())->nadirVector();
    TVector<double> ideal = static_cast<BBOBBiObj*>(func.get())->idealVector();
    double bestHyp = static_cast<BBOBBiObj*>(func.get())->bestHypervolumeValue();

    form->defineProblem(prob);

    init->TP_defineSetSize(popSize);

    alg->defineReferenceSetSize(nDirs);
    alg->defineScalarisingFunction(sFunc);
    alg->TP_defineMaxSolutions(maxSurrogateSize);
    alg->defineOptimizationSearchQuality(0);

    tend->defineBudget(budget);

    LogManagerSPtr log = tend->log();
    log->defineLogAll(true);

    coco_archive_t* arc = coco_archive("bbob-biobj", functionID, nInputs, instanceID);

    TVector<double> Icoco;
    TVector<double> IcocoComp;
    bool firstIteration = true;
    int nSolsArchive = 0;

    using std::chrono::steady_clock;
    steady_clock::time_point start = steady_clock::now();

    while(tend->remainingBudget() > 0) {
        tend->evaluate();

        TVector<IMappingSPtr> evaluatedMappings = alg->evaluatedMappings();

        if(firstIteration) {
            TVector<IMappingSPtr> currentEvalMappings;
            for(size_t i=0; i<evaluatedMappings.size(); i++) {
                IMappingSPtr sol = evaluatedMappings[i];
                TString sName = "s" + std::to_string(i);
                TVector<double> z = sol->doubleObjectiveVec();
                coco_archive_add_solution(arc, z[0], z[1], sName.data());

                currentEvalMappings.push_back(sol);
                nSolsArchive = coco_archive_get_number_of_solutions(arc);
                if(nSolsArchive == 2) {
                    double minDist = minDistRI(currentEvalMappings, ideal, nadir);
                    Icoco.push_back(minDist);
                    IcocoComp.push_back(minDist);
                } else {
                    double hpv = -coco_archive_get_hypervolume(arc);
                    Icoco.push_back(hpv);
                    ISet* archive = base->setWithTag(Tigon::TNonDominatedArchive);
                    hpv = -hypervolumeTigon(archive, ideal, nadir);
                    IcocoComp.push_back(hpv);
                }
            }

            firstIteration = false;
        } else {
            IMappingSPtr sol = evaluatedMappings.back();
            TVector<double> z = sol->doubleObjectiveVec();
            TString sName = "s" + std::to_string(evaluatedMappings.size());
            int nondominated = coco_archive_add_solution(arc, z[0], z[1], sName.data());
            if(nondominated) {
                nSolsArchive = coco_archive_get_number_of_solutions(arc);
            }
            if(nSolsArchive == 2) {
                double minDist = minDistRI(evaluatedMappings, ideal, nadir);
                Icoco.push_back(minDist);
                IcocoComp.push_back(minDist);
            } else {
                double hpv = -coco_archive_get_hypervolume(arc);
                Icoco.push_back(hpv);
                ISet* archive = base->setWithTag(Tigon::TNonDominatedArchive);
                hpv = -hypervolumeTigon(archive, ideal, nadir);
                IcocoComp.push_back(hpv);
            }
        }

        cout << "Reamining budget: " << tend->remainingBudget() << endl;
        cout << " Direction vector: \t";
        dispVector(alg->dirVec(), " ", "\n");
        cout << " Number of IMappings: "
             << alg->evaluatedMappings().size() + 1 << endl;

        tend->incrementIteration();
    }
    steady_clock::time_point end   = steady_clock::now();

    ISet* pop = new ISet(alg->evaluatedMappings());

#ifdef LOG_POPULATION
    // Save population and evaluations
    JsonObject customFields;
    customFields["Random Seed"] = seed;

    JsonArray nadirArray = {nadir[0], nadir[1]};
    customFields["Nadir Vector"] = nadirArray;

    JsonArray idealArray = {ideal[0], ideal[1]};
    customFields["Ideal Vector"] = idealArray;

    customFields["Best Hypervolume"] = bestHyp;

    customFields["TimeMS"] = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    JsonArray hypervolumeArray;
    for(double hpv : Icoco) {
        hypervolumeArray.push_back(hpv);
    }
    customFields["ICoco"] = hypervolumeArray;

    JsonArray hypervolumeArrayComp;
    for(double hpv : IcocoComp) {
        hypervolumeArrayComp.push_back(hpv);
    }
    customFields["ICocoComp"] = hypervolumeArrayComp;

    log->logPopulation(pop, "Main Population", customFields);

    TString filePath = "outputs/ParEGO/D" + std::to_string(nInputs) + "/runD" + std::to_string(nInputs) + "F" + std::to_string(functionID);
    std::filesystem::create_directories(filePath);

    JsonDocument jdoc;
    jdoc.setObject(log->populationsLog());
    TString fileName;
    fileName = filePath + "/population_ParEGO_BBOB-BIOBJ_F"
               + std::to_string(functionID) + "_D"
               + std::to_string(nInputs) + "_I"
               + std::to_string(instanceID) + "_S"
               + std::to_string(seed) + ".json";
    Tigon::writeToFile(TString(jdoc.toJson().c_str()), fileName);

    JsonDocument jdoc2;
    jdoc2.setObject(log->evaluationsLog());
    TString fileName2;
    fileName2 = filePath + "/evaluations_ParEGO_BBOB-BIOBJ_F"
                + std::to_string(functionID) + "_D"
                + std::to_string(nInputs) + "_I"
                + std::to_string(instanceID) + "_S"
                + std::to_string(seed) + ".json";
    Tigon::writeToFile(TString(jdoc2.toJson().c_str()), fileName2);
#endif

    coco_archive_free(arc);

    delete tend;
    delete alg;
    delete evaluator;
    delete init;
    delete form;
    delete base;
}

void bmk_MParEGO_BBOBBIOBJ(int budget, int seed, int functionID, int nInputs, int instanceID, 
                           Tigon::ScalarisationType sFunc=Tigon::WeightedChebyshevAugmented, 
                           int numberSamples=200)
{
    TRAND.defineSeed(seed);
    cout.precision(5);

    int   nDirs    = 11;
    int   popSize  = 11 * nInputs - 1;
    // int   maxSurrogateSize = popSize;
    int   maxSurrogateSize = 100;

    PSetBase*     base      = new PSetBase();
    IFormulation* form      = new IFormulation(base);
    SobolInit*    init      = new SobolInit(form);
    Evaluator*    evaluator = new Evaluator(init);
    MultiParEGO*  alg       = new MultiParEGO(evaluator);
    ITermination* tend      = new ITermination(alg);

    base->defineKeepArchive(true);

    bool useObserver = false;
    IFunctionSPtr func = BiObjProblem(functionID, nInputs, instanceID,
                                      useObserver, "MultiParEGO");
    ProblemSPtr prob(new Problem);
    prob->appendFunction(func);
    prob->processProblemDefinition();

    TVector<double> nadir = static_cast<BBOBBiObj*>(func.get())->nadirVector();
    TVector<double> ideal = static_cast<BBOBBiObj*>(func.get())->idealVector();
    double bestHyp = static_cast<BBOBBiObj*>(func.get())->bestHypervolumeValue();

    form->defineProblem(prob);

    init->TP_defineSetSize(popSize);

    alg->defineReferenceSetSize(nDirs);
    alg->defineScalarisingFunction(sFunc);
    alg->defineNumberSamples(numberSamples);
    alg->TP_defineMaxSolutions(maxSurrogateSize);
    alg->defineOptimizationSearchQuality(0);

    tend->defineBudget(budget);

    LogManagerSPtr log = tend->log();
    log->defineLogAll(true);

    coco_archive_t* arc = coco_archive("bbob-biobj", functionID, nInputs, instanceID);

    TVector<double> Icoco;
    TVector<double> IcocoComp;
    bool firstIteration = true;
    int nSolsArchive = 0;

    using std::chrono::steady_clock;
    steady_clock::time_point start = steady_clock::now();

    while(tend->remainingBudget() > 0) {
        tend->evaluate();

        TVector<IMappingSPtr> evaluatedMappings = alg->evaluatedMappings();

        if(firstIteration) {
            TVector<IMappingSPtr> currentEvalMappings;
            for(size_t i=0; i<evaluatedMappings.size(); i++) {
                IMappingSPtr sol = evaluatedMappings[i];
                TString sName = "s" + std::to_string(i);
                TVector<double> z = sol->doubleObjectiveVec();
                coco_archive_add_solution(arc, z[0], z[1], sName.data());

                currentEvalMappings.push_back(sol);
                nSolsArchive = coco_archive_get_number_of_solutions(arc);
                if(nSolsArchive == 2) {
                    double minDist = minDistRI(currentEvalMappings, ideal, nadir);
                    Icoco.push_back(minDist);
                    IcocoComp.push_back(minDist);
                } else {
                    double hpv = -coco_archive_get_hypervolume(arc);
                    Icoco.push_back(hpv);
                    ISet* archive = base->setWithTag(Tigon::TNonDominatedArchive);
                    hpv = -hypervolumeTigon(archive, ideal, nadir);
                    IcocoComp.push_back(hpv);
                }
            }

            firstIteration = false;
        } else {
            IMappingSPtr sol = evaluatedMappings.back();
            TVector<double> z = sol->doubleObjectiveVec();
            TString sName = "s" + std::to_string(evaluatedMappings.size());
            int nondominated = coco_archive_add_solution(arc, z[0], z[1], sName.data());
            if(nondominated) {
                nSolsArchive = coco_archive_get_number_of_solutions(arc);
            }
            if(nSolsArchive == 2) {
                double minDist = minDistRI(evaluatedMappings, ideal, nadir);
                Icoco.push_back(minDist);
                IcocoComp.push_back(minDist);
            } else {
                double hpv = -coco_archive_get_hypervolume(arc);
                Icoco.push_back(hpv);
                ISet* archive = base->setWithTag(Tigon::TNonDominatedArchive);
                hpv = -hypervolumeTigon(archive, ideal, nadir);
                IcocoComp.push_back(hpv);
            }
        }

        cout << "Reamining budget: " << tend->remainingBudget() << endl;
        cout << " Direction vector: \t";
        dispVector(alg->dirVec(), " ", "\n");
        cout << " Number of IMappings: "
             << alg->evaluatedMappings().size() + 1 << endl;

        tend->incrementIteration();
    }
    steady_clock::time_point end   = steady_clock::now();

    ISet* pop = new ISet(alg->evaluatedMappings());

#ifdef LOG_POPULATION
    // Save population and evaluations
    JsonObject customFields;
    customFields["Random Seed"] = seed;

    JsonArray nadirArray = {nadir[0], nadir[1]};
    customFields["Nadir Vector"] = nadirArray;

    JsonArray idealArray = {ideal[0], ideal[1]};
    customFields["Ideal Vector"] = idealArray;

    customFields["Best Hypervolume"] = bestHyp;

    customFields["TimeMS"] = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    JsonArray hypervolumeArray;
    for(double hpv : Icoco) {
        hypervolumeArray.push_back(hpv);
    }
    customFields["ICoco"] = hypervolumeArray;

    JsonArray hypervolumeArrayComp;
    for(double hpv : IcocoComp) {
        hypervolumeArrayComp.push_back(hpv);
    }
    customFields["ICocoComp"] = hypervolumeArrayComp;

    log->logPopulation(pop, "Main Population", customFields);

    TString filePath = "outputs/MParEGO/D" + std::to_string(nInputs) + "/runD" + std::to_string(nInputs) + "F" + std::to_string(functionID);
    std::filesystem::create_directories(filePath);

    JsonDocument jdoc;
    jdoc.setObject(log->populationsLog());
    TString fileName;
    fileName = filePath + "/population_MParEGO_BBOB-BIOBJ_F"
               + std::to_string(functionID) + "_D"
               + std::to_string(nInputs) + "_I"
               + std::to_string(instanceID) + "_S"
               + std::to_string(seed) + ".json";
    Tigon::writeToFile(TString(jdoc.toJson().c_str()), fileName);

    JsonDocument jdoc2;
    jdoc2.setObject(log->evaluationsLog());
    TString fileName2;
    fileName2 = filePath + "/evaluations_MParEGO_BBOB-BIOBJ_F"
                + std::to_string(functionID) + "_D"
                + std::to_string(nInputs) + "_I"
                + std::to_string(instanceID) + "_S"
                + std::to_string(seed) + ".json";
    Tigon::writeToFile(TString(jdoc2.toJson().c_str()), fileName2);
#endif

    coco_archive_free(arc);

    delete tend;
    delete alg;
    delete evaluator;
    delete init;
    delete form;
    delete base;
}

int main(int argc, char* argv[]) 
{
    using std::chrono::steady_clock;

    // defauls
    std::string algoName = "Random"; // Random | Sobol | ParEGO | MParEGO
    int functionID = 54;
    int nInputs = 3;
    int budget = 1000;
    int numberSamples = 200;
    Tigon::ScalarisationType sFunc=Tigon::WeightedChebyshevAugmented;

    // simple CLI parsing:
    // argv[1] = algorithm name
    // argv[2] = budget
    // argv[3] = functionID
    // argv[4] = nInputs
    // argv[5] = numberSamples (only used by MParEGO)

    if(argc >= 2) algoName = argv[1];
    if(argc >= 3) budget = atoi(argv[2]);
    if(argc >= 4) functionID = atoi(argv[3]);
    if(argc >= 5) nInputs = atoi(argv[4]);
    if(argc >= 6) numberSamples = atoi(argv[5]);

    auto print_usage = [&](const char* prog){
        std::cout << "Usage: " << prog << " <algorithm> [budget] [functionID] [nInputs] [numberSamples]\n"
                  << "  algorithm: Random | Sobol | ParEGO | MParEGO\n"
                  << "  budget default: " << budget << "\n"
                  << "  functionID default: " << functionID << "\n"
                  << "  nInputs default: " << nInputs << "\n"
                  << "  numberSamples default (MParEGO only): " << numberSamples << "\n";
    };

    if(algoName != "Random" && algoName != "Sobol" && algoName != "ParEGO" && algoName != "MParEGO") {
        print_usage(argv[0]);
        return -1;
    }

    steady_clock::time_point start = steady_clock::now();

    // int instanceID=1;
    for(int instanceID=1; instanceID<=15; instanceID++) {
        if(algoName == "Random") {
            bmk_Random_BBOBBIOBJ(budget, instanceID, functionID, nInputs, instanceID);
        } else if(algoName == "Sobol") {
            bmk_Sobol_BBOBBIOBJ(budget, instanceID, functionID, nInputs, instanceID);
        } else if(algoName == "ParEGO") {
            bmk_ParEGO_BBOBBIOBJ(budget, instanceID, functionID, nInputs, instanceID, sFunc);
        } else if(algoName == "MParEGO") {
            bmk_MParEGO_BBOBBIOBJ(budget, instanceID, functionID, nInputs, instanceID, sFunc, numberSamples);
        } else {
            print_usage(argv[0]);
            return -1;
        }
    }

    steady_clock::time_point end   = steady_clock::now();

    std::cout << "Time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms" <<std::endl;

    return 0;
}
