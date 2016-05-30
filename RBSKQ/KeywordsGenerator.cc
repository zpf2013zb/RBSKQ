//
//  KeywordsGenerator.cpp
//  rksk_query
//
//  Created by 秦旭 on 14-7-1.
//  Copyright (c) 2014年 ZJU. All rights reserved.
//

#include "KeywordsGenerator.h"
#include <bitset>
#include <vector>
#include <time.h>
//#define MAX_KEYWORDS 64
#define MAX_KEYWORDS 64

//std::vector<unsigned long long> KeywordsGenerator::getKeywords(std::size_t totalNumberOfNodes, std::size_t avgKeywords, std::vector<unsigned long long>& keys)
//{
//	srand((unsigned)time(NULL));
//	std::default_random_engine generator;
//	std::normal_distribution<double> distribution(avgKeywords, 1);
//
//	std::vector<unsigned long long> tmp;
//	//tmp.clear();
//	tmp.reserve(totalNumberOfNodes);
//	std::bitset<MAX_KEYWORDS> keywordsSet;
//
//	while (totalNumberOfNodes--) {
//		keywordsSet.reset();
//		int numberOfKeywords = 0;
//		while (numberOfKeywords < 1) {
//			numberOfKeywords = distribution(generator);
//		}
//		//int percent = random(100)%1000;
//		int percent = rand() % 1000;
//		// precent used for what?
//		if (percent < 50) {
//
//			while (keywordsSet.count() < numberOfKeywords)
//			{
//				// int index = random()%64;
//				int index = rand() % 64;
//				keywordsSet.set(index);
//
//			}
//
//		}
//		else {
//
//			while (keywordsSet.count() < numberOfKeywords)
//			{
//				int index;
//				if (keywordsSet.count() < 6) {
//
//					// index = random()%8;
//					index = rand() % 8;
//				}
//				else {
//					//index = random() % 20;
//					index = rand() % 20;
//				}
//
//
//				keywordsSet.set(index + 22);
//
//			}
//
//		}
//
//		tmp.push_back(keywordsSet.to_ullong());
//
//	}
//	keys.assign(tmp.begin(), tmp.end());
//	return tmp;
//}


std::vector<unsigned long long> KeywordsGenerator::getKeywords(std::size_t totalNumberOfNodes, std::size_t avgKeywords, std::vector<unsigned long long>& keys)
{
	srand((unsigned)time(NULL));
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(avgKeywords, 1);

	std::vector<unsigned long long> tmp;
	//tmp.clear();
	tmp.reserve(totalNumberOfNodes);
	std::bitset<MAX_KEYWORDS> keywordsSet;

	while (totalNumberOfNodes--) {
		keywordsSet.reset();
		int numberOfKeywords = 0;
		while (numberOfKeywords < 1) {
			numberOfKeywords = distribution(generator);
		}
		int percent = rand() % 1000;
		// precent used for what?
		if (percent < 100) {

			while (keywordsSet.count() < numberOfKeywords)
			{
				// int index = random()%64;
				int index = rand() % 64;
				keywordsSet.set(index);

			}

		}
		else {

			while (keywordsSet.count() < numberOfKeywords)
			{
				int index;
				if (keywordsSet.count() < 6) {

					// index = random()%8;
					index = rand() % 16;
				}
				else {
					//index = random() % 20;
					index = rand() % 32;
				}


				keywordsSet.set(index);

			}

		}

		tmp.push_back(keywordsSet.to_ullong());

	}
	keys.assign(tmp.begin(), tmp.end());
	return tmp;
}


//int KeywordsGenerator::getNumOfPoints(int aveNPt)
//{
//	srand((unsigned)time(NULL));
//	std::default_random_engine generator;
//	std::normal_distribution<double> distribution(aveNPt, 1);
//
//	int numberOfKeywords = 0;
//	while (numberOfKeywords < 1) {
//		numberOfKeywords = distribution(generator);
//	}
//	return numberOfKeywords;
//}

int KeywordsGenerator::getNumOfPoints(int aveNPt)
{
	srand((unsigned)time(NULL));
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(aveNPt, 1);

	int numberOfKeywords = 0;
	while (numberOfKeywords < 1) {
		numberOfKeywords = distribution(generator);
	}
	return numberOfKeywords;
}


//std::vector<unsigned long long> KeywordsGenerator::getConstantKeywords(std::size_t totalNumberOfNodes, std::size_t number)
//{
//    std::vector<unsigned long long> keys;
//    keys.reserve(totalNumberOfNodes);
//    std::bitset<MAX_KEYWORDS> keywords_set;
//	srand((unsigned)time(NULL));
//    while (totalNumberOfNodes --) {
//        
//        keywords_set.reset();
//        
//        while(keywords_set.count() < number )
//        {
//            int index;
//            if (keywords_set.count() < 6) {
//                //index = random()%8;
//				index = rand()%8;
//            }
//			else {
//				//index = random() % 20;
//				index = rand() % 20;
//			}
//
//            keywords_set.set(index+22);
//            
//        }//Avoid generate vct#=0
//        
//        keys.push_back(keywords_set.to_ullong());
//    }
//    
//    return keys;
//    
//}

std::vector<unsigned long long> KeywordsGenerator::getConstantKeywords(std::size_t totalNumberOfNodes, std::size_t number)
{
	std::vector<unsigned long long> keys;
	keys.reserve(totalNumberOfNodes);
	std::bitset<MAX_KEYWORDS> keywords_set;
	srand((unsigned)time(NULL));
	while (totalNumberOfNodes--) {

		keywords_set.reset();
		int percent = rand() % 1000;
		if (percent < 100) {

			while (keywords_set.count() < number)
			{
				// int index = random()%64;
				int index = rand() % 64;
				keywords_set.set(index);

			}

		}
		else {

			while (keywords_set.count() < number)
			{
				int index;
				if (keywords_set.count() < 6) {
					// index = random()%8;
					index = rand() % 16;
				}
				else {
					//index = random() % 20;
					index = rand() % 32;
				}
				keywords_set.set(index);

			}

		}

		keys.push_back(keywords_set.to_ullong());
	}
	return keys;
}


KeywordsGenerator::KeywordsGenerator()
{
}

KeywordsGenerator::~KeywordsGenerator()
{
}


