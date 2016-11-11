#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <limits>
#include <cfloat>
using namespace std;

/*

CS 180 Programming Project
Spring 2014
Patrick Tan
Disscussion 1B

*/

/*void priA(double ** printArray, int col, int row) //for debugging purposes
{
	for(int i = 0; i < col; i++)
	{
		for(int j = 0; j< row; j++)
		{
			cout << printArray[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void priB(double * printArray, int n) //for debugging purposes
{
	for(int i = 0; i < n; i++)
	{
		cout << printArray[i] << " ";
		cout << endl;
	}
	cout << endl;
}*/

int main(int argc, char *argv[])
{
	//Get input file
	ifstream input(argv[1]);
	//ifstream input;
	//input.open("test.txt");
	if(!input.is_open())
	{
		cout << "Error Opening";
		return 0;
	}

	//Take lecture length and number of topics
	int lectureLength;
	int numTopics;
	double temp;
	input >> lectureLength;
	input >> numTopics;

	//Create topic array and fill with topic lengths
	//O(n)
	double* topicLengths = new double[numTopics];
	for(int i = 0; i < numTopics; i++)
	{
		input >> temp;
		topicLengths[i] = temp;
	}

	//create empty two dimensional array of lectures for memoization
	//O(n^2)
	double **lectures = new double*[numTopics];
	for(int i = 0; i < numTopics; i++)
	{
		lectures[i] = new double[numTopics];
	}

	/*
	This is a multiway choice problem.

	OPT(n) = dis(i,n) + OPT(i-1)
	for each lecture, we want to find the optimal dissatisfaction with optimal solution above

	We can solve it in a similar way as the matrix example in class.

	We will create an array of the dissatisfaction.
	Then, we will find the array for a cutoff for topics to get the least dissatisfaction using 
	OPT(j) = min(1<=i<=j) (dis(i,n) + OPT(i-1))
		segment pi to pj is used in an optimum solution if it is part of the the minimum obtained using index i
	Finally, we can output the best set of lectures using the both arrays
	*/

	/*
	lectures array will hold the amount of dissatisfaction for a given set of topics in a single lecture
	rows will be the starting topic and columns will be the ending topic
	i.e. column 2 row 2 will mean the time leftover if only topic three is in the lecture
		 column 2 row 3 will be the time left over if topics three and four are both done in the lecture
	the upper right triangle will be empty because a the topics must be done sequentially
	if the number goes into the negatives, this indicates this topic cannot possibly be done in this lecture (not enough time)

	O(n^2)
	*/

	//Run through all sets of lectures to set number of excess minutes
	for(int i = 0; i < numTopics; i++)
	{
		for(int j = 0; j < numTopics; j++)
		{
			if(i < j)
				lectures[i][j] = 0;
			else if (i == j)
			{
				lectures[i][j] = lectureLength - topicLengths[i];
			}
			else
			{
				lectures[i][j] = lectures[i - 1][j] - topicLengths[i];
			}
		}
	}
	
	//priA(lectures, numTopics, numTopics);

	//Given number of excess minutes, set dissatisfaction
	double temporary;
	for(int i = 0; i < numTopics; i++)
	{
		for(int j = 0; j < numTopics; j++)
		{
			if(i == numTopics - 1 && lectures[i][j] >= 0) //if it's the last lecture, no dissatisfaction
				lectures[i][j] = 0;
			else if(lectures[i][j] <= 5 && lectures[i][j] >= 0) //if less than five minutes excess, no dissatisfaction
				lectures[i][j] = 0;
			else if(lectures[i][j] > 5)
			{
				temporary = pow((double)lectures[i][j] - 5, (double)4);
				if(temporary >= DBL_MAX)
					temporary = DBL_MAX;
				lectures[i][j] = temporary; //dissatisfaction for extra time
			}
			else if(lectures[i][j] < 0) //-1 signifies over time limit
				lectures[i][j] = -1;
		}
	}
		
	//priA(lectures, numTopics, numTopics);

	/*
	Now, we use the recurrence to find the minimum dissatisfaction of given subsets.

	We use the recurrence for the multiway solution
	OPT(j) = min(1<=i<=j) (lectures[i][j] + optimal[i-1])
		segment pi to pj is used in an optimum solution if it is part of the the minimum obtained using index i

	We make two arrays, optimal and dissatisfaction
	Optimal holds the optimal set of topics in a subset and dissatisfaction holds total dissatisfactions for a given subset.
	We will use both for the output.

	O(n^2)
	*/

	double* optimal = new double[numTopics];
	double* dissatisfaction = new double[numTopics];
	dissatisfaction[0] = lectures[0][0];
	optimal[0] = 0;
	double tempMin = DBL_MAX;
	double minDis = DBL_MAX;
	for(int j = 1; j < numTopics; j++)
	{
		if(lectures[j][0] != -1)
		{
			minDis = lectures[j][0];
			tempMin = 0;
		}
		for(int i = 1; i <= j; i++)
		{
			//lectures 2 through 4 would be lectures[4][2]
			if(lectures[j][i] != -1)
			{
				if(lectures[j][i] + dissatisfaction[i - 1] < minDis)
				{
					minDis = lectures[j][i] + dissatisfaction[i - 1];
					tempMin = i;
				}
			}
		}
		dissatisfaction[j] = minDis;
		optimal[j] = tempMin;

		tempMin = DBL_MAX;
		minDis = DBL_MAX;
	}

	//priB(dissatisfaction, numTopics);
	//priB(optimal, numTopics);

	/*
	We use optimal and dissatisfaction to output the results.
	We have to traverse optimal in reverse order, so to output in the correct order, 
		we create some new array output where we put the results. We then go through 
		this array and output everything in the correct order.

	O(n)
	*/

	int j = numTopics - 1;
	int numLectures = 0;
	double * output = new double[numTopics];

	while(j >= 0)
	{
		output[numLectures] = optimal[j] + 1;
		j = optimal[j] - 1;
		numLectures++;
	}

	//priB(output, numTopics);
	cout << fixed << setprecision(0);
	cout << numLectures << endl;

	numLectures--;

	while(numLectures > 0)
	{
		cout << output[numLectures] << "-" << output[numLectures - 1] - 1 << endl;
		numLectures--;
	}

	cout << output[0] << "-" << numTopics << endl;

	cout << dissatisfaction[numTopics - 1] << endl;

	//deleting allocated arrays
	for(int i = 0; i < numTopics; i++)
	{
		delete [] lectures[i];
	}
	delete [] lectures;
	delete [] topicLengths;
	delete [] optimal;
	delete [] dissatisfaction;
	delete [] output;
}