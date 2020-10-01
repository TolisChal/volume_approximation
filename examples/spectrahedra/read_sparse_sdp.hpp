

#ifndef READ_SPARSE_SDP_HPP
#define READ_SPARSE_SDP_HPP

typedef std::string::iterator string_it3;
typedef std::vector<double> listVector3;

char consumeSymbol3(string_it3 &at, string_it3 &end) {
    while (at != end) {
        if (*at != ' ' && *at != '\t') {
            char c = *at;
            at++;
            return c;
        }

        at++;
    }

    return '\0';
}


bool isCommentLine3(std::string &line) {
    string_it3 at = line.begin();
    string_it3 end = line.end();

    char c = consumeSymbol3(at, end);

    return c == '"' || c == '*';
}


int fetchNumber3(std::string &string) {
    std::stringstream stream(string);
    int num;
    stream >> num;
    return num;
}

std::vector<double> readVector3(std::string &string) {
    std::stringstream stream(string);
    std::vector<double> vector;
    double value;

    while (stream >> value) {
        vector.push_back(value);
    }

    return vector;
}


/// Reads an SDPA format file
/// \param[in] is An open stram pointing to the file
/// \param[out] matrices the matrices A0, A1, A2, ..., An
/// \param[out] objectiveFunction The objective function of the sdp
template <typename MT, typename LMII, typename VT>
void loadSparseSDPAFormatFile3(std::istream &is, LMII &lmi, VT &objectiveFunction) {

    //std::cout << "hello1" << "\n";
    std::string line;
    std::string::size_type sz;

    std::getline(is, line, '\n');

    //skip comments
    while (isCommentLine3(line)) {
        std::getline(is, line, '\n');
    }

    //read variables number
    int variablesNum = fetchNumber3(line);
    //std::cout << "variablesNum: " << variablesNum <<std::endl;

    if (std::getline(is, line, '\n').eof())
        throw std::runtime_error("Unexpected end of file");

    //read number of blocks
    int blocksNum = fetchNumber3(line);

    if (std::getline(is, line, '\n').eof())
        throw std::runtime_error("Unexpected end of file");

    //read block structure vector
    listVector3 blockSizes = readVector3(line);

    if (blockSizes.size() != blocksNum)
        throw std::runtime_error("Wrong number of blocks");

    if (std::getline(is, line, '\n').eof())
        throw std::runtime_error("Unexpected end of file");

    //read objective function
    listVector3 constantVector = readVector3(line);

    //std::cout << "hello2" << "\n";
    while (constantVector.size() < variablesNum) {
        if (std::getline(is, line, '\n').eof())
            throw std::runtime_error("Unexpected end of file");

        listVector3 t = readVector3(line);
        constantVector.insert(std::end(constantVector), std::begin(t), std::end(t));
    }
    //std::cout << "hello3" << "\n";

    std::vector<MT> matrices = std::vector<MT>(variablesNum + 1);
    int matrixDim = 0;
    for (auto x : blockSizes)
        matrixDim += std::abs((int) x);

    //std::cout << "hello4" << "\n";
    for (int i=0 ; i<matrices.size() ; ++i)
        matrices[i].setZero(matrixDim, matrixDim);
    // read constraint matrices
    // entries are of the form
    // <matno> <blkno> <i> <j> <entry>

    //std::cout << "hello5" << "\n";
    while (!std::getline(is, line, '\n').eof()) {
        listVector3 t = readVector3(line);
            //std::cout << line << "\n";
        int blockOffset = 0;
        for (int i=1; i<t[1] ; ++i) blockOffset += std::abs(blockSizes[i-1]);

        int i = t[2] + blockOffset-1;
        int j = t[3] + blockOffset-1;

        matrices[t[0]](i,j) = t[4];
//            std::cout << i << " " << j << "\n";
        // matrix is symmetric
        // only upper triangular is provided
        // fill lower triangular
        if (i!=j) matrices[t[0]](j,i) = t[4];
    }
    //std::cout << "hello6" << "\n";
    for (int atMatrix=1 ; atMatrix<matrices.size() ; atMatrix++) {
        //the LMI in SDPA format is >0, I want it <0
        //F0 has - before it in SDPA format, the rest have +
        matrices[atMatrix] *= -1;
    }
    //std::cout << "hello7" << "\n";
    // return lmi and objective function
    objectiveFunction.setZero(variablesNum);
    int at = 0;

    //std::cout << "matrices.size(): " << matrices.size() <<std::endl;
    for (auto value : constantVector)
        objectiveFunction(at++) = value;
    lmi = LMII(matrices);
}

#endif
