//
// Created by Hao Shi on 8/21/17.
//

#include "../include/MDCas2s.h"

using namespace std;
using namespace tensor_hao;

MDCas2s::MDCas2s():L(0), Nup(0), Ndn(0), MDSize(0), MDUpSize(0), MDDnSize(0), NupMax(0), NdnMax(0), logw(0.0), wfCasIsCalculated(false) {}

MDCas2s::MDCas2s(const string &filename) { read(filename); }

MDCas2s::MDCas2s(const MDCas2s &x) { copy_deep(x); }

MDCas2s::MDCas2s(MDCas2s &&x) { move_deep(x); }

MDCas2s::~MDCas2s() { }

MDCas2s &MDCas2s::operator=(const MDCas2s &x) { copy_deep(x); return *this; }

MDCas2s &MDCas2s::operator=(MDCas2s &&x) { move_deep(x); return *this; }

size_t MDCas2s::getL() const { return L; }

size_t MDCas2s::getNup() const { return Nup; }

size_t MDCas2s::getNdn() const { return Ndn; }

size_t MDCas2s::getMDSize() const { return MDSize; }

size_t MDCas2s::getMDUpSize() const { return MDUpSize; }

size_t MDCas2s::getMDDnSize() const { return MDDnSize; }

size_t MDCas2s::getNupMax() const { return NupMax; }

size_t MDCas2s::getNdnMax() const { return NdnMax; }

const complex<double> &MDCas2s::getLogw() const { return logw; }

const vector<complex<double>> &MDCas2s::getCoe() const { return coe; }

const vector<size_t> &MDCas2s::getCoeLinkUp() const { return coeLinkUp; }

const vector<size_t> &MDCas2s::getCoeLinkDn() const { return coeLinkDn; }

const vector<vector<size_t>> &MDCas2s::getOccupancyUp() const { return occupancyUp; }

const vector<vector<size_t>> &MDCas2s::getOccupancyDn() const { return occupancyDn; }

const TensorHao<complex<double>, 2> &MDCas2s::getWfUp() const { return wfUp; }

const TensorHao<complex<double>, 2> &MDCas2s::getWfDn() const { return wfDn; }

const TensorHao<complex<double>, 2> &MDCas2s::getRUp() const { return RUp; }

const TensorHao<complex<double>, 2> &MDCas2s::getRDn() const { return RDn; }

const vector<int> &MDCas2s::getIsParentUp() const { return isParentUp; }

const vector<int> &MDCas2s::getIsParentDn() const { return isParentDn; }

const vector<size_t> &MDCas2s::getTreeUp() const { return treeUp; }

const vector<size_t> &MDCas2s::getTreeDn() const { return treeDn; }

const vector<vector<size_t>> &MDCas2s::getDiffOrbitalUp() const { return diffOrbitalUp; }

const vector<vector<size_t>> &MDCas2s::getDiffOrbitalDn() const { return diffOrbitalDn; }

const vector<size_t> &MDCas2s::getParentIndexUp() const { return parentIndexUp; }

const vector<size_t> &MDCas2s::getParentIndexDn() const { return parentIndexDn; }

const TensorHao<complex<double>, 2> &MDCas2s::generateWfCasUp() { calculateWfCas(); return wfCasUp; }

const TensorHao<complex<double>, 2> &MDCas2s::generateWfCasDn() { calculateWfCas(); return wfCasDn; }

void MDCas2s::stabilize()
{
    calculateWfCas();

    wfUp = wfCasUp;
    wfDn = wfCasDn;

    BL_NAME(QRMatrix)(wfUp, RUp);
    BL_NAME(QRMatrix)(wfDn, RDn);

    scaleRMatrix();
}

complex<double> MDCas2s::normalize()
{
    stabilize();
    complex<double> logwTemp(logw);
    logw=0.0;
    return logwTemp;
}

void MDCas2s::read(const std::string &filename)
{
    ifstream file;
    file.open(filename, ios::in);
    if ( ! file.is_open() ) { cout<<"Error opening file in File!!! "<<filename<<endl; exit(1); }

    readFile(L, file);      readFile(Nup, file);      readFile(Ndn, file);
    readFile(MDSize, file); readFile(MDUpSize, file); readFile(MDDnSize, file); readFile(NupMax, file); readFile(NdnMax, file);
    readFile(logw, file);
    coe.resize(MDSize); readFile(coe.size(), coe.data(), file);
    coeLinkUp.resize(MDSize); readFile(coeLinkUp.size(), coeLinkUp.data(), file);
    coeLinkDn.resize(MDSize); readFile(coeLinkDn.size(), coeLinkDn.data(), file);
    occupancyUp.resize(MDUpSize);
    for(size_t i = 0; i < MDUpSize; ++i)
    {
        occupancyUp[i].resize(Nup);
        readFile(Nup, occupancyUp[i].data(), file);
    }
    occupancyDn.resize(MDDnSize);
    for(size_t i = 0; i < MDDnSize; ++i)
    {
        occupancyDn[i].resize(Ndn);
        readFile(Ndn, occupancyDn[i].data(), file);
    }
    wfUp.resize(L, NupMax); readFile(wfUp.size(), wfUp.data(), file);
    wfDn.resize(L, NdnMax); readFile(wfDn.size(), wfDn.data(), file);
    RUp.resize(NupMax, NupMax); readFile( RUp.size(), RUp.data(), file );
    RDn.resize(NdnMax, NdnMax); readFile( RDn.size(), RDn.data(), file );

    file.close();

    //scaleRMatrix need to use wfCasIsCalculated
    wfCasIsCalculated=false;

    normCoe();
    scaleRMatrix();
    setUpParentLinkStructure();
    setParentIndex();

    printParentNumber();
}

void MDCas2s::write(const std::string &filename) const
{
    ofstream file;
    file.open(filename, ios::out|ios::trunc);
    if ( ! file.is_open() ) { cout << "Error opening file in File!!! "<<filename<<endl; exit(1); }

    writeFile(L, file);      writeFile(Nup, file);      writeFile(Ndn, file);
    writeFile(MDSize, file); writeFile(MDUpSize, file); writeFile(MDDnSize, file); writeFile(NupMax, file); writeFile(NdnMax, file);
    writeFile(logw, file);

    vector< complex<double> > coeOrigin(MDSize);
    vector<int> signUp, signDn;
    vector< vector<size_t> > occupancyUpOrigin, occupancyDnOrigin;

    signUp = sortOrbitalsInOccupancy(occupancyUpOrigin, occupancyUp);
    signDn = sortOrbitalsInOccupancy(occupancyDnOrigin, occupancyDn);
    for(size_t i = 0; i < MDSize ; ++i)
    {
        coeOrigin[i] = signUp[ coeLinkUp[i] ] * signDn[ coeLinkDn[i] ] * 1.0 * coe[i];
    }

    writeFile(coeOrigin.size(), coeOrigin.data(), file);
    writeFile(coeLinkUp.size(), coeLinkUp.data(), file);
    writeFile(coeLinkDn.size(), coeLinkDn.data(), file);
    for(size_t i = 0; i < MDUpSize; ++i) writeFile(Nup, occupancyUpOrigin[i].data(), file);
    for(size_t i = 0; i < MDDnSize; ++i) writeFile(Ndn, occupancyDnOrigin[i].data(), file);
    writeFile(wfUp.size(), wfUp.data(), file);
    writeFile(wfDn.size(), wfDn.data(), file);
    writeFile(RUp.size(), RUp.data(), file);
    writeFile(RDn.size(), RDn.data(), file);

    file.close();
}

void MDCas2s::printParentNumber()
{
    size_t parentNumber;

    cout<<"\nTotal number of determinant is "<<MDSize<<"."<<endl;

    parentNumber=0; for(size_t i = 0; i <MDUpSize; ++i) parentNumber+= isParentUp[i];
    cout<<"Number of determinant for spin up is "<<MDUpSize<<"."<<endl;
    cout<<"Number of parent for spin up is "<<parentNumber<<"."<<endl;

    parentNumber=0; for(size_t i = 0; i <MDDnSize; ++i) parentNumber+= isParentDn[i];
    cout<<"Number of determinant for spin dn is "<<MDDnSize<<"."<<endl;
    cout<<"Number of parent for spin dn is "<<parentNumber<<".\n"<<endl;
}

double MDCas2s::getMemory() const
{
    double mem(0.0);
    mem += 8.0*8+16.0;
    mem += 16.0*coe.size();
    mem += 8.0*coeLinkUp.size()+8.0*coeLinkDn.size();
    mem += 8.0*occupancyUp.size()*Nup + 8.0*occupancyDn.size()*Ndn;
    mem += wfUp.getMemory() + wfDn.getMemory();
    mem += RUp.getMemory() + RDn.getMemory();
    mem += 4.0*isParentUp.size() + 4.0*isParentDn.size();
    mem += 8.0*treeUp.size() + 8.0*treeDn.size();
    for(size_t i = 0; i < MDUpSize; ++i) mem += 8.0*diffOrbitalUp[i].size();
    for(size_t i = 0; i < MDDnSize; ++i) mem += 8.0*diffOrbitalDn[i].size();
    mem += 8.0*parentIndexUp.size() + 8.0*parentIndexDn.size();
    mem += 1.0 + wfCasUp.getMemory() + wfCasDn.getMemory();
    return mem;
}

#ifdef MPI_HAO
void MPIBcast(MDCas2s &buffer, int root, MPI_Comm const &comm)
{
    MPIBcast( buffer.L, root, comm );
    MPIBcast( buffer.Nup, root, comm );
    MPIBcast( buffer.Ndn, root, comm );
    MPIBcast( buffer.MDSize, root, comm );
    MPIBcast( buffer.MDUpSize, root, comm );
    MPIBcast( buffer.MDDnSize, root, comm );
    MPIBcast( buffer.NupMax, root, comm );
    MPIBcast( buffer.NdnMax, root, comm );
    MPIBcast( buffer.logw, root, comm );
    MPIBcast( buffer.coe, root, comm );
    MPIBcast( buffer.coeLinkUp, root, comm );
    MPIBcast( buffer.coeLinkDn, root, comm );

    buffer.occupancyUp.resize(buffer.MDUpSize);
    for(size_t i = 0; i < buffer.MDUpSize; ++i)
    {
        MPIBcast( buffer.occupancyUp[i], root, comm );
    }
    buffer.occupancyDn.resize(buffer.MDDnSize);
    for(size_t i = 0; i < buffer.MDDnSize; ++i)
    {
        MPIBcast( buffer.occupancyDn[i], root, comm );
    }
    MPIBcast( buffer.wfUp, root, comm );
    MPIBcast( buffer.wfDn, root, comm );
    MPIBcast( buffer.RUp, root, comm );
    MPIBcast( buffer.RDn, root, comm );

    MPIBcast( buffer.isParentUp, root, comm );
    MPIBcast( buffer.isParentDn, root, comm );
    MPIBcast( buffer.treeUp, root, comm );
    MPIBcast( buffer.treeDn, root, comm );

    buffer.diffOrbitalUp.resize(buffer.MDUpSize);
    for(size_t i = 0; i < buffer.MDUpSize; ++i)
    {
        MPIBcast( buffer.diffOrbitalUp[i], root, comm );
    }
    buffer.diffOrbitalDn.resize(buffer.MDDnSize);
    for(size_t i = 0; i < buffer.MDDnSize; ++i)
    {
        MPIBcast( buffer.diffOrbitalDn[i], root, comm );
    }
    MPIBcast( buffer.parentIndexUp, root, comm );
    MPIBcast( buffer.parentIndexDn, root, comm );
    MPIBcast( buffer.wfCasIsCalculated, root, comm );
    MPIBcast( buffer.wfCasUp, root, comm );
    MPIBcast( buffer.wfCasDn, root, comm );
}
#endif

void MDCas2s::copy_deep(const MDCas2s &x)
{
    L = x.L;
    Nup = x.Nup;
    Ndn = x.Ndn;
    MDSize = x.MDSize;
    MDUpSize = x.MDUpSize;
    MDDnSize = x.MDDnSize;
    NupMax = x.NupMax;
    NdnMax = x.NdnMax;
    logw = x.logw;
    coe = x.coe;
    coeLinkUp = x.coeLinkUp;
    coeLinkDn = x.coeLinkDn;
    occupancyUp = x.occupancyUp;
    occupancyDn = x.occupancyDn;
    wfUp = x.wfUp;
    wfDn = x.wfDn;
    RUp = x.RUp;
    RDn = x.RDn;
    isParentUp = x.isParentUp;
    isParentDn = x.isParentDn;
    treeUp = x.treeUp;
    treeDn = x.treeDn;
    diffOrbitalUp = x.diffOrbitalUp;
    diffOrbitalDn = x.diffOrbitalDn;
    parentIndexUp = x.parentIndexUp;
    parentIndexDn = x.parentIndexDn;
    wfCasIsCalculated = x.wfCasIsCalculated;
    wfCasUp = x.wfCasUp;
    wfCasDn = x.wfCasDn;
}

void MDCas2s::move_deep(MDCas2s &x)
{
    L = x.L;
    Nup = x.Nup;
    Ndn = x.Ndn;
    MDSize = x.MDSize;
    MDUpSize = x.MDUpSize;
    MDDnSize = x.MDDnSize;
    NupMax = x.NupMax;
    NdnMax = x.NdnMax;
    logw = x.logw;
    coe = move( x.coe );
    coeLinkUp = move(x.coeLinkUp);
    coeLinkDn = move(x.coeLinkDn);
    occupancyUp = move( x.occupancyUp );
    occupancyDn = move( x.occupancyDn );
    wfUp = move( x.wfUp );
    wfDn = move( x.wfDn );
    RUp = move( x.RUp );
    RDn = move( x.RDn );
    isParentUp = move( x.isParentUp );
    isParentDn = move( x.isParentDn );
    treeUp = move( x.treeUp );
    treeDn = move( x.treeDn );
    diffOrbitalUp = move( x.diffOrbitalUp );
    diffOrbitalDn = move( x.diffOrbitalDn );
    parentIndexUp = move( x.parentIndexUp );
    parentIndexDn = move( x.parentIndexDn );
    wfCasIsCalculated = x.wfCasIsCalculated;
    wfCasUp = move( x.wfCasUp );
    wfCasDn = move( x.wfCasDn );
}

void MDCas2s::calculateWfCas()
{
    if( wfCasIsCalculated ) return;

    wfCasUp.resize(L, NupMax); BL_NAME(gmm)(wfUp, RUp, wfCasUp);
    wfCasDn.resize(L, NdnMax); BL_NAME(gmm)(wfDn, RDn, wfCasDn);

    wfCasIsCalculated = true;
}

void MDCas2s::normCoe()
{
    double norm(0.0);
    for(size_t i = 0; i < MDSize; ++i) norm += abs( coe[i] ) * abs( coe[i] );
    norm = sqrt(norm);
    for(size_t i = 0; i < MDSize; ++i) coe[i] /= norm;
    logw += log( norm );
}

void MDCas2s::scaleRMatrix()
{
    double maxValue;

    maxValue = -1.0;
    for(size_t i = 0; i < NupMax; ++i)
    {
        if( maxValue < abs( RUp(i,i) ) ) maxValue = abs( RUp(i,i) );
    }
    RUp /= complex<double>(maxValue);
    if( wfCasIsCalculated ) wfCasUp /= complex<double>(maxValue);
    logw += Nup * log(maxValue);

    maxValue = -1.0;
    for(size_t i = 0; i < NdnMax; ++i)
    {
        if( maxValue < abs( RDn(i,i) ) ) maxValue = abs( RDn(i,i) );
    }
    RDn /= complex<double>(maxValue);
    if( wfCasIsCalculated ) wfCasDn /= complex<double>(maxValue);
    logw += Ndn * log(maxValue);
}

void MDCas2s::setUpParentLinkStructure()
{
    int isOrder;
    vector<int> signUp, signDn;

    if(MDUpSize==0) { cout<<"Error!!! Number of determinant for spin up can not be zero!"<<endl; exit(1); }
    if(MDDnSize==0) { cout<<"Error!!! Number of determinant for spin dn can not be zero!"<<endl; exit(1); }

    isOrder = checkOccupancyOrder(occupancyUp);
    if( !isOrder ) { cout<<"Error!!! The multi-determiant orbital for spin up is not in order!"<<endl; exit(1); }
    buildParentTree(isParentUp, treeUp, occupancyUp);
    signUp = permuteOrbitalsInOccupancy(isParentUp, treeUp, occupancyUp);
    diffOrbitalUp = calculateDiffOrbitalsInOccupancy(isParentUp, treeUp, occupancyUp);

    isOrder = checkOccupancyOrder(occupancyDn);
    if( !isOrder ) { cout<<"Error!!! The multi-determiant orbital for spin dn is not in order!"<<endl; exit(1); }
    buildParentTree(isParentDn, treeDn, occupancyDn);
    signDn = permuteOrbitalsInOccupancy(isParentDn, treeDn, occupancyDn);
    diffOrbitalDn = calculateDiffOrbitalsInOccupancy(isParentDn, treeDn, occupancyDn);

    for(size_t i = 0; i < MDSize; ++i)
    {
        coe[i] *= ( signUp[ coeLinkUp[i] ] * signDn[ coeLinkDn[i] ] * 1.0 );
    }
}

void MDCas2s::setParentIndex()
{
    size_t currentParent, currentIndex;

    currentParent = 0; parentIndexUp.resize(MDUpSize);
    for(size_t i = 0; i < MDUpSize; ++i)
    {
        currentIndex = treeUp[i];
        parentIndexUp[currentIndex] = currentParent;
        if( isParentUp[currentIndex] == 1 ) currentParent = currentIndex;
    }

    currentParent = 0; parentIndexDn.resize(MDDnSize);
    for(size_t i = 0; i < MDDnSize; ++i)
    {
        currentIndex = treeDn[i];
        parentIndexDn[currentIndex] = currentParent;
        if( isParentDn[currentIndex] == 1 ) currentParent = currentIndex;
    }
}