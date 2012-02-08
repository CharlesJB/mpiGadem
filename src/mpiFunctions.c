
void broadcast2DimensionsDoubleArray(double* doublePtr, int numberOfRows, int elementsPerRow, int source)
void broadcast1DimensionIntArray(int* intPtr, int numberOfElements, source);
void broadcast1DimensionDoubleArray(double* doublePtr, int numberOfElements, source);
void broadcast1DimensionCharArray(char* charPtr, int numberOfElements, source);

broadcastOnce(int maxSeqLen, 
              int numEM, 
              int startPWMfound, 
              int minminSites, 
              double maxpFactor, 
              int numSeq, 
              int numSeqEM, 
              char* Iseq;
              double* bfreq0, 
              double** posWeight, 
              int weightType, 
              double pvalueCutoff, 
              int* emSeqLen, 
              int populationSize)
{
    int* intPtr = NULL;
    double* doublePtr = NULL;
    int source = 0;

    broadcast1DimensionIntArray(&maxSeqLen, 1, source);
    broadcast1DimensionIntArray(&numEM, 1, source);
    broadcast1DimensionIntArray(&startPWMfound, 1, source);
    broadcast1DimensionIntArray(&minminSites, 1, source);
    broadcast1DimensionDoubleArray(&maxpFactor, 1, source);
    broadcast1DimensionIntArray(&numSeq, 1, source);
    broadcast1DimensionIntArray(&numSeqEM, 1, source);
    broadcast1DimensionCharArray(Iseq, numSeq+1, source);
    broadcast1DimensionDoubleArray(bfreq0, 4, source);
    broadcast2DimensionsDoubleArray(posWeight, numSeq, maxSeqLen, source); 
    broadcast1DimensionIntArray(&weightType, 1, source);
    /* pValueCutoff */
    broadcast1DimensionIntArray(emSeqLen, numSeqEM, source);
    
       
}
broadcastEveryCycle(double** pwm, 
                  int pwmLen, 
                  char* pwmConsensus, 
                  int scoreCutoff, 
                  char* sdyad[ii], 
                  int populationSize)

    int* intPtr = NULL;
    double* doublePtr = NULL;
    char* charPtr = NULL;
    int numberOfElements = 0;
    int source = 0;
    int i = 0;
    double* doubleIterator = NULL;

    broadcast3DimensionsDoubleArray(pwm, populationSize, MAX_PWM_LENGTH, 4, source)
    broadcast1DimensionIntArray(&pwmLen, populationSize, source);
    
    /* pwmConsensus */
    broadcast1DimensionIntArray(&scoreCutoff, /* ??? */, source);
    /* sdyad */       
    
    
}

void broadcast1DimensionIntArray(int* intPtr, int numberOfElements, int source)
{
    MPI_Bcast(intPtr, numberOfElements, MPI_INT, source, MPI_COMM_WORLD);
}

void broadcast1DimensionDoubleArray(double* doublePtr, int numberOfElements, int source)
{
    MPI_Bcast(doublePtr, numberOfElements, MPI_DOUBLE, source, MPI_COMM_WORLD);
}

void broadcast1DimensionCharArray(char* charPtr, int numberOfElements, int source)
{
    MPI_Bcast(doublePtr, numberOfElements, MPI_CHAR, source, MPI_COMM_WORLD);
}

void broadcast2DimensionsDoubleArray(double* doublePtr, int numberOfRows, int elementsPerRow, int source)
{
    int i = 0; 
    int size[2];
    double* doubleIterator = doublePtr;
    
    size[0] = numberOfRows;
    size[1] = elementsPerRow;    

    broadcast1DimensionIntArray(size, 2, source);

    for (i = 0; i < numberOfRows; i++)
    {
        broadcast1DimensionDoubleArray(doubleIterator, elementsPerRow, source);
        doubleIterator = doubleIterator + elementsPerRow;
    }
}

void broadcast2DimensionsCharArray(char* charPtr, int numberOfRows, int elementsPerRow, int source)
{
    int i = 0; 
    int size[2];
    char* charIterator = charPtr;
    
    size[0] = numberOfRows;
    size[1] = elementsPerRow;    

    broadcast1DimensionIntArray(size, 2, source);

    for (i = 0; i < numberOfRows; i++)
    {
        broadcast1DimensionCharArray(charIterator, elementsPerRow, source); 
        charIterator = charIterator + elementsPerRow;
    }
}

void broadcast3DimensionsDoubleArray(double* doublePtr, int numberOfArrays, int numberOfRowsPerArray, int elementsPerRow, int source)
{
    int i - 0;
    double* doubleIterator = NULL;

    broadcast1DimensionIntArray(&numberOfArrays, 1, source);
    doubleIterator = doublePtr;

    for (i = 0; i < populationSize; i++)

void broadcast3DimensionsDoubleArray(double* doublePtr, int numberOfArrays, int numberOfRowsPerArray, int elementsPerRow, int source)
{
    int i - 0;
    double* doubleIterator = NULL;

    broadcast1DimensionIntArray(&numberOfArrays, 1, source);
    doubleIterator = doublePtr;

    for (i = 0; i < populationSize; i++)
    {
        broadcast2DimensionsDoubleArray(doubleIterator, numberOfRowsPerArray, elementsPerRow, source);
        doubleIterator = doubleIterator + (numberOfRowsPerArray*elementsPerRow);
    }
}

