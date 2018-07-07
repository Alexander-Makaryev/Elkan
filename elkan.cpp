#include <math.h>

//the most expensive operation, we try to reduce the number of calls
double distance(int p, double* x, double* y) {
	double temp;
	double d = 0;
	for(int i = 0; i < p; i++) {
		temp = x[i]-y[i];
		d += temp * temp;
	}
	return sqrt(d);
}

//the names of variables are similar to the article by Charles Elkan
//The numbers in the comments are steps from the article
void elkan(int p, int n, double** x, int k, double** cStart, double threshold, int maxIterationsCount,
double** c, int* pointsCenters, int* iterationsCount, bool* isEmptyCluster) {
	double temp, temp2;
	int* cPointsCount = new int[k];
	bool* r = new bool[n];
	double* u = new double[n];
	double* s = new double[k];
	bool* centerChanged = new bool[k];
	double* centersDelta = new double[k];
	int* pointsMaxDisIndex = new int[k];

	double** l = new double*[n];
	for(int i = 0; i < n; i++) {
		l[i] = new double[k];
	}

	double** m = new double*[k];
	for(int i = 0; i < k; i++) {
		m[i] = new double[p];
	}
	
	double** cd = new double*[k];
	for(int i = 0; i < k; i++) {
		cd[i] = new double[k];
	}

	*isEmptyCluster = false;

	for(int i = 0; i < k; i++) {
		for(int j = 0; j < p; j++) {
			c[i][j] = cStart[i][j];
		}
	}
	
	for(int i = 0; i < k; i++) {
		for(int j = i+1; j < k; j++) {
			temp = distance(p, c[i], c[j]);
			cd[i][j] = temp;
			cd[j][i] = temp;
		}
	}

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < k; j++) {
			l[i][j] = 0;
		}
	}

	for(int i = 0; i < n; i++) {
		pointsCenters[i] = 0;
		u[i] = distance(p, x[i], c[0]);
		l[i][0] = u[i];
		for(int j = 1; j < k; j++) {
			if(u[i] > cd[pointsCenters[i]][j]/2.0) {
				temp = distance(p, x[i], c[j]);
				l[i][j] = temp;
				if(temp < u[i]) {
					pointsCenters[i] = j;
					u[i] = temp;
				}
			}
		}
	}


	bool isCentersChanged = true;
	int currentIterationNumber = 0;
	while((currentIterationNumber < maxIterationsCount) && isCentersChanged) {
		currentIterationNumber++;

		//1
		bool isValue;
		for(int i = 0; i < k; i++) {
			isValue = false;
			for(int j = 0; j < k; j++) {
				if(i != j) {
					if(isValue) {
						if(cd[i][j] < temp) {
							temp = cd[i][j];
						}
					} else {
						temp = cd[i][j];
						isValue = true;
					}
				}
			}
			s[i] = temp/2.0;
		}

		//2&3
		for(int i = 0; i < n; i++) {
			if(u[i] > s[pointsCenters[i]]) {
				for(int j = 0; j < k; j++) {
					if((pointsCenters[i] != j) &&
						(u[i] > l[i][j]) &&
						(u[i] > cd[j][pointsCenters[i]]/2.0)
					) {
						if(r[i]) {
							temp = distance(p, x[i], c[pointsCenters[i]]);
							l[i][pointsCenters[i]] = temp;
                            u[i] = temp;
							r[i] = false;
						} else {
							temp = u[i];
						}
					
                        if((temp > l[i][j]) || (temp > cd[pointsCenters[i]][j]/2.0)) {
                            temp2 = distance(p, x[i], c[j]);
                            l[i][j] = temp2;
                            if(temp2 < temp) {
                                pointsCenters[i] = j;
                                u[i] = temp2;
                            }
                        }
                    }
				}
			}
		}

		//4
		for(int i = 0; i < k; i++) {
			cPointsCount[i] = 0;
			for(int j = 0; j < p; j++) {
				m[i][j] = 0;
			}
		}
		for(int i = 0; i < n; i++) {
			cPointsCount[pointsCenters[i]]++;
			for(int j = 0; j < p; j++) {
				m[pointsCenters[i]][j] += x[i][j];
			}
		}
		int emptyClustersCount = 0;
		for(int i = 0; i < k; i++) {
			if (cPointsCount[i] > 0) {
				for(int j = 0; j < p; j++) {
					m[i][j] = m[i][j] / (double)cPointsCount[i];
				}
			} else {
				emptyClustersCount++;
			}
		}

		if(emptyClustersCount > 0) {
			*isEmptyCluster = true;

			int minMaxIndex = 0;
			for(int i = 0; i < emptyClustersCount; i++) {
				pointsMaxDisIndex[i] = i;
				if (u[i] < u[minMaxIndex]) {
					minMaxIndex = i;
				}
			}
			for(int i = emptyClustersCount; i < n; i++) {
				if(u[i] > u[pointsMaxDisIndex[minMaxIndex]]) {
					pointsMaxDisIndex[minMaxIndex] = i;
					minMaxIndex = 0;
					for(int j = 1; j < emptyClustersCount; j++) {
						if (u[pointsMaxDisIndex[j]] < u[pointsMaxDisIndex[minMaxIndex]]) {
							minMaxIndex = j;
						}
					}
				}
			}
			int tempInt;
			for(int i = 0; i < emptyClustersCount; i++) {
				for(int j = i+1; j < emptyClustersCount; j++) {
					if (u[pointsMaxDisIndex[i]] < u[pointsMaxDisIndex[j]]) {
						tempInt = pointsMaxDisIndex[i];
						pointsMaxDisIndex[i] = pointsMaxDisIndex[j];
						pointsMaxDisIndex[j] = tempInt;
					}
				}
			}
			
			tempInt = 0;
			for(int i = 0; i < k; i++) {
				if(cPointsCount[i] == 0) {
					for(int j = 0; j < p; j++) {
						m[i][j] = x[pointsMaxDisIndex[tempInt]][j];
					}
					tempInt++;
				}
			}
		}

		//threshold check
		temp = 0;
		for(int i = 0; i < k; i++) {
            centersDelta[i] = distance(p, c[i], m[i]);
			temp += centersDelta[i];
        }
		if(temp * temp < threshold) {
			isCentersChanged = false;
		}

		//detect changed centers
		for(int i = 0; i < k; i++) {
			centerChanged[i] = false;
			for(int j = 0; j < p; j++) {
				if(c[i][j] != m[i][j]) {
					centerChanged[i] = true;
					break;
				}
			}
		}
		//only if there is a next iteration
		if(isCentersChanged) {
			//compute the distances between centers(from step 1)
			for(int i = 0; i < k; i++) {
				if(centerChanged[i]) {
					for(int j = 0; j < k; j++) {
						if((i != j) &&
							!((j < i) && centerChanged[j])
							) {
							temp = distance(p, m[i], m[j]);
							cd[i][j] = temp;
							cd[j][i] = temp;
						}
					}
				}
			}

			//5
			for(int i = 0; i < n; i++) {
				for(int j = 0; j < k; j++) {
					temp = l[i][j] - centersDelta[j];
					l[i][j] = temp > 0 ? temp : 0;
				}
			}
			//6
			for(int i = 0; i < n; i++) {
				u[i] = u[i] + centersDelta[pointsCenters[i]];
				r[i] = true;
			}
		}
		//7
		for(int i = 0; i < k; i++) {
			if(centerChanged[i]) {
				for(int j = 0; j < p; j++) {
					c[i][j] = m[i][j];
				}
			}
		}
	}

	delete[] cPointsCount;
	delete[] r;
	delete[] u;
	delete[] s;
	delete[] centerChanged;
	delete[] centersDelta;
	delete[] pointsMaxDisIndex;
	for(int i = 0; i < n; i++) {
		delete[] l[i];
	}
	delete[] l;
	for(int i = 0; i < k; i++) {
		delete[] m[i];
	}
	delete[] m;
	for(int i = 0; i < k; i++) {
		delete[] cd[i];
	}
	delete[] cd;

	*iterationsCount = currentIterationNumber;
}

