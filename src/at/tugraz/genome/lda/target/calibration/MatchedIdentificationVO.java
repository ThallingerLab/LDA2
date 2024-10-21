/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */

package at.tugraz.genome.lda.target.calibration;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.commons.math3.exception.OutOfRangeException;

import org.apache.commons.math3.util.Pair;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class MatchedIdentificationVO
{
	private ArrayList<IdentificationVO> originals_;
	private int maxOriginals_;
	private ArrayList<IdentificationVO> matches_;
	private int maxMatches_;
	private RecalibrationRegression regression_;
	private double predictionThreshold_;
	private String matchingAlgo_;
	private Pair<IdentificationVO,IdentificationVO> highestConfidencePair_;
	private ArrayList<Pair<IdentificationVO,IdentificationVO>> acceptedMatches_ = new ArrayList<Pair<IdentificationVO,IdentificationVO>>();
	private int confidence_ = 0;
//	ArrayList<Pair<IdentificationVO,IdentificationVO>> elutionOrderAssignments_ = new ArrayList<Pair<IdentificationVO,IdentificationVO>>();
	
	protected MatchedIdentificationVO(ArrayList<IdentificationVO> originals, int maxOriginals, ArrayList<IdentificationVO> matches, int maxMatches, 
			RecalibrationRegression regression, double predictionThreshold, String matchingAlgo)
	{
		this.originals_ = originals;
		this.maxOriginals_ = maxOriginals;
		this.matches_ = matches;
		this.maxMatches_ = maxMatches;
		this.regression_ = regression;
		this.predictionThreshold_ = predictionThreshold;
		this.matchingAlgo_ = matchingAlgo;
		findHighestConfidencePairings();
//		findElutionOrderAssignments();
	}
	
	/**
	 * The algorithm finds pairs as follows:
	 * Take the identification with the highest average area for both original and new conditions => highest confidence pair (confidence = 0)
	 * (confidence 0 is always accepted for standards, 1 is required for biological data).
	 * 
	 * If the matching algorithm is CalibrationFileChooserPanel.MATCHING_SETTING_INCREASED_RELIABILITY or
	 * the standards regression is null or the retention times are outside the range of the regression and the matching algorithm is CalibrationFileChooserPanel.MATCHING_SETTING_DEFAULT: 
	 * 1. If this identification is found in all given results files
	 * 2. If this identification has MSn evidence in at least one file
	 * 3. If the number of MSn peaks is inconsistent (1 and more than 1), we risk false matches. => confidence remains 0
	 * 4. If there are more than one MSn peaks for both conditions: if the elution order matches the intensity order
	 * then this identification (only the one with the highest intensity) is matched with confidence = 1.
	 * 
	 * Else (CalibrationFileChooserPanel.MATCHING_SETTING_INCREASED_NUMBER):
	 * 1. Take identifications with available MSn evidence, which are found in all given result files
	 * 2. Find the closest identification to the predicted retention time for each
	 * 3. Take preferentially hits with highest total area (of original and new measurements combined)
	 * 4. Check if the retention time of the matches are consistent for both conditions
	 */
	private void findHighestConfidencePairings()
	{
		Collections.sort(this.originals_, obtainAreaComparator());
		Collections.sort(this.matches_, obtainAreaComparator());
		IdentificationVO highestIntensityOriginal = originals_.get(0);
		IdentificationVO highestIntensityMatch = matches_.get(0);
		this.highestConfidencePair_ = new Pair<IdentificationVO,IdentificationVO>(highestIntensityOriginal, highestIntensityMatch);
		
		if (matchingAlgo_.equals(CalibrationFileChooserPanel.MATCHING_SETTING_ALL_MATCHES)) //all data recieves the minimum confidence required for being accepted
		{
			this.confidence_ = 1;
		}
		
		if ((isWithinRTRange() || matchingAlgo_.equals(CalibrationFileChooserPanel.MATCHING_SETTING_INCREASED_NUMBER)) &&
				!matchingAlgo_.equals(CalibrationFileChooserPanel.MATCHING_SETTING_INCREASED_RELIABILITY))
		{
			this.acceptedMatches_ = findAcceptedMatches(sortPredictedMatchesByTotalArea());
			if (!acceptedMatches_.isEmpty())
			{
				this.confidence_ = 2;
			}
		}
		else
		{
			if (isAlwaysFound(highestIntensityOriginal, maxOriginals_) && 
					isAlwaysFound(highestIntensityMatch, maxMatches_) &&
					isIntensityElutionOrderCorrelated())
			{
				this.confidence_ = 1;
			}
		}
	}
	
	/**
	 * Finds accepted matches relying on considerations based on the area (matches with a higher total area will be chosen).
	 * If the retention time behaviour in the original and new measurement is consistent, the found matches are returned (otherwise an empty ArrayList).
	 * @param matchesArea
	 * @return
	 */
	private ArrayList<Pair<IdentificationVO,IdentificationVO>> findAcceptedMatches(PriorityQueue<Pair<IdentificationVO,IdentificationVO>> matchesArea)
	{
		Set<IdentificationVO> originalSet = ConcurrentHashMap.newKeySet();
		Set<IdentificationVO> newSet = ConcurrentHashMap.newKeySet();
		PriorityQueue<Pair<IdentificationVO,IdentificationVO>> matchesRTKey = new PriorityQueue<Pair<IdentificationVO,IdentificationVO>>(obtainRTKeyComparator());
		PriorityQueue<Pair<IdentificationVO,IdentificationVO>> matchesRTValue = new PriorityQueue<Pair<IdentificationVO,IdentificationVO>>(obtainRTValueComparator());
		ArrayList<Pair<IdentificationVO,IdentificationVO>> matchesAccepted = new ArrayList<Pair<IdentificationVO,IdentificationVO>>();
		while (!matchesArea.isEmpty())
		{
			Pair<IdentificationVO,IdentificationVO> match = matchesArea.poll();
			if (!originalSet.contains(match.getKey()) && !newSet.contains(match.getValue()))
			{
				matchesRTKey.add(match);
				matchesRTValue.add(match);
				originalSet.add(match.getKey());
				newSet.add(match.getValue());
			}
		}
		boolean isRTConsistent = true;
		while (!matchesRTKey.isEmpty() && !matchesRTValue.isEmpty())
		{
			Pair<IdentificationVO,IdentificationVO> matchRTKey = matchesRTKey.poll();
			matchesAccepted.add(matchRTKey);
			if (!matchRTKey.equals(matchesRTValue.poll()))
			{
				isRTConsistent = false;
			}
		}
		if (!isRTConsistent)
		{
			matchesAccepted.clear();
		}
		return matchesAccepted;
	}
	
	/**
	 * After obtaining hypothetical matches based on the regression obtained from standards, 
	 * this method sorts these matches by the combined area of the original and new measurements in a priority queue.
	 * @return
	 */
	private PriorityQueue<Pair<IdentificationVO,IdentificationVO>> sortPredictedMatchesByTotalArea()
	{
		Set<IdentificationVO> originalSet = ConcurrentHashMap.newKeySet();
		Set<IdentificationVO> newSet = ConcurrentHashMap.newKeySet();
		PriorityQueue<Pair<IdentificationVO,IdentificationVO>> matchesPrediction = new PriorityQueue<Pair<IdentificationVO,IdentificationVO>>(obtainPredictionComparator());
		PriorityQueue<Pair<IdentificationVO,IdentificationVO>> matchesArea = new PriorityQueue<Pair<IdentificationVO,IdentificationVO>>(obtainAreasComparator());
		
		for (IdentificationVO originalVO : originals_)
		{
			if (originalVO.isMSnAvailable() && isAlwaysFound(originalVO,maxOriginals_))
			{
				originalSet.add(originalVO);
				for (IdentificationVO newVO : matches_)
				{
					if (newVO.isMSnAvailable() && isAlwaysFound(newVO,maxMatches_))
					{
						if (matchingAlgo_.equals(CalibrationFileChooserPanel.MATCHING_SETTING_INCREASED_NUMBER) || isWithinPredictionThreshold(originalVO, newVO))
						{
							newSet.add(newVO);
							matchesPrediction.add(new Pair<IdentificationVO,IdentificationVO>(originalVO, newVO));
						}
					}
				}
			}
		}
		while ((!originalSet.isEmpty() || !newSet.isEmpty()) && !matchesPrediction.isEmpty())
		{
			Pair<IdentificationVO,IdentificationVO> match = matchesPrediction.poll();
			originalSet.remove(match.getKey());
			newSet.remove(match.getValue());
			matchesArea.add(match);
		}
		return matchesArea;
	}
	
	/**
	 * Checks whether the retention time difference between identifications measured on the two different conditions
	 * is within the defined threshold.
	 * @param originalVO
	 * @param newVO
	 * @return
	 */
	private boolean isWithinPredictionThreshold(IdentificationVO originalVO, IdentificationVO newVO)
	{
		double predictedRT = regression_.getTargetRT(originalVO.getAverageRT());
		double predictionDiff = Math.abs(predictedRT - newVO.getAverageRT());
		if (predictionDiff < predictionThreshold_)
		{
			return true;
		}
		return false;
	}
	
	/**
	 * Checks for each identification if the retention time is within the range of the regression. 
	 * @return
	 */
	private boolean isWithinRTRange()
	{
		if (regression_ != null)
		{
			for (IdentificationVO vo : originals_)
			{
				try
				{
					regression_.getTargetRT(vo.getAverageRT());
				}
				catch (OutOfRangeException ex) 
				{
					return false;
				}
				catch (Exception ex)
				{
					ex.printStackTrace();
				}
			}
			return true;
		}
		return false;
	}
	
	/**
	 * This method assumes that originals_ and matches_ are ordered in descending order by their average area.
	 * If there is only one MSn peak for both, the method returns true. If the number of MSn peaks is inconsistent (1 and more than 1), the method returns false.
	 * If there are at least two MSn peaks, it checks whether the elution order matches the intensity order.
	 * @return
	 */
	private boolean isIntensityElutionOrderCorrelated()
	{
		if (getNumberMSn(originals_) == 1 && getNumberMSn(matches_) == 1 &&
				areBothMSn(originals_.get(0), matches_.get(0)))
		{
			return true; //if the highest intensity peak is the one and only MSn peak, we trust it.
		}
		else if (getNumberMSn(originals_) == 1 || getNumberMSn(matches_) == 1)
		{
			return false; //if the number of MSn peaks is inconsistent (1 and more than 1), we risk false matches.
		}
		else if (originals_.size() > 1 && matches_.size() > 1 &&
				     isCorrelatedTopDown())
		{
			return true; //if both the highest and second highest intensity peaks are MSn hits and the elution order is preserved, we trust it.
		}
		else
		{
			return false;
		}
	}
	
	/**
	 * This method assumes that originals_ and matches_ are ordered in descending order by their average area.
	 * It checks whether the elution order matches the intensity order for the minimum number of MSn peaks of originals_ and matches_.
	 * Starting with the third highest peaks, if not both are MSn but the criteria were fulfilled thus far, the method returns true.
	 * @return
	 */
	private boolean isCorrelatedTopDown()
	{
		boolean isCorrelated = false;
		for (int i=0; i<Math.min(originals_.size(), matches_.size())-1; i++)
		{
			if (areBothMSn(originals_.get(i), matches_.get(i)) &&
			    areBothMSn(originals_.get(i+1), matches_.get(i+1)))
			{
				if ((originals_.get(i).getAverageRT() > originals_.get(i+1).getAverageRT() && matches_.get(i).getAverageRT() > matches_.get(i+1).getAverageRT()) 
							 ||
					  (originals_.get(i).getAverageRT() < originals_.get(i+1).getAverageRT() && matches_.get(i).getAverageRT() < matches_.get(i+1).getAverageRT()))
				{
					isCorrelated = true;
				}
				else
				{
					return false; //when a mismatch is detected, we return immediately.
				}
			}
			else if (i == 0) //the first two identifications must be MSn.
			{
				return false;
			}
			else
			{
				return isCorrelated;
			}
		}
		return isCorrelated;
	}
	
	private int getNumberMSn(ArrayList<IdentificationVO> vos)
	{
		int count = 0;
		for (IdentificationVO vo : vos)
		{
			if (vo.isMSnAvailable())
			{
				count ++;
			}
		}
		return count;
	}
	
	private boolean isAlwaysFound(IdentificationVO vo, int max)
	{
		return vo.getIdentificationCount() == max;
	}
	
	private boolean areBothMSn(IdentificationVO vo1, IdentificationVO vo2)
	{
		return vo1.isMSnAvailable() && vo2.isMSnAvailable();
	}
	
	private Comparator<IdentificationVO> obtainAreaComparator()
	{
		return new Comparator<IdentificationVO>() 
		{
			@Override
			public int compare(IdentificationVO o1, IdentificationVO o2)
			{
				return - o1.getAverageArea().compareTo(o2.getAverageArea());
			}
		};
	}
	
	private Comparator<Pair<IdentificationVO,IdentificationVO>> obtainAreasComparator()
	{
		return new Comparator<Pair<IdentificationVO,IdentificationVO>>() 
		{
			@Override
			public int compare(Pair<IdentificationVO,IdentificationVO> o1, Pair<IdentificationVO,IdentificationVO> o2)
			{
				return - new Double(o1.getKey().getAverageArea() + o1.getValue().getAverageArea())
						.compareTo(o2.getKey().getAverageArea() + o2.getValue().getAverageArea());
			}
		};
	}
	
	private Comparator<Pair<IdentificationVO,IdentificationVO>> obtainRTKeyComparator()
	{
		return new Comparator<Pair<IdentificationVO,IdentificationVO>>() 
		{
			@Override
			public int compare(Pair<IdentificationVO,IdentificationVO> o1, Pair<IdentificationVO,IdentificationVO> o2)
			{
				return o1.getKey().getAverageRT().compareTo(o2.getKey().getAverageRT());
			}
		};
	}
	
	private Comparator<Pair<IdentificationVO,IdentificationVO>> obtainRTValueComparator()
	{
		return new Comparator<Pair<IdentificationVO,IdentificationVO>>() 
		{
			@Override
			public int compare(Pair<IdentificationVO,IdentificationVO> o1, Pair<IdentificationVO,IdentificationVO> o2)
			{
				return o1.getValue().getAverageRT().compareTo(o2.getValue().getAverageRT());
			}
		};
	}
	
	private Comparator<Pair<IdentificationVO,IdentificationVO>> obtainPredictionComparator()
	{
		return new Comparator<Pair<IdentificationVO,IdentificationVO>>() 
		{
			@Override
			public int compare(Pair<IdentificationVO,IdentificationVO> o1, Pair<IdentificationVO,IdentificationVO> o2)
			{
				if (!matchingAlgo_.equals(CalibrationFileChooserPanel.MATCHING_SETTING_INCREASED_NUMBER))
				{
					double predictedRT1 = regression_.getTargetRT(o1.getKey().getAverageRT());
					double predictedRT2 = regression_.getTargetRT(o2.getKey().getAverageRT());
					
					return new Double(Math.abs(predictedRT1 - o1.getValue().getAverageRT()))
							.compareTo(Math.abs(predictedRT2 - o2.getValue().getAverageRT()));
				}
				else //the setting for increased number of identifications also works without a standards-curve for predicted values
				{
					return new Double(o1.getValue().getAverageRT())
							.compareTo(o2.getValue().getAverageRT());
				}
			}
		};
	}
	
	public ArrayList<Pair<IdentificationVO,IdentificationVO>> getAcceptedMatches()
	{
		return acceptedMatches_;
	}

	public Pair<IdentificationVO,IdentificationVO> getHighestConfidencePair()
	{
		return highestConfidencePair_;
	}

	public int getConfidence()
	{
		return confidence_;
	}
	
	/**
	 * Has to be executed after the highest confidence pair is found.
	 * Identifications are sorted based on retention time. Then they are added to an ArrayList
	 * first descending and then ascending starting from the highest confidence pair.
	 */
//	private void findElutionOrderAssignments()
//	{
//		//TODO: make sure the sorting is done in the correct order!
//		Collections.sort(originals_, obtainRTComparator());
//		Collections.sort(matches_, obtainRTComparator());
//		int indexOriginal = originals_.indexOf(highestConfidencePair_.getKey());
//		int indexMatch = matches_.indexOf(highestConfidencePair_.getValue());
//		
//		int j = indexMatch;
//		for (int i = indexOriginal-1; i>=0; i--)
//		{
//			j--;
//			if (j < 0) continue;
//			elutionOrderAssignments_.add(new Pair<IdentificationVO,IdentificationVO>(originals_.get(i),matches_.get(j)));
//		}
//		j = indexMatch;
//		for (int i = indexOriginal+1; i<originals_.size(); i++)
//		{
//			j++;
//			if (j >= matches_.size()) continue;
//			elutionOrderAssignments_.add(new Pair<IdentificationVO,IdentificationVO>(originals_.get(i),matches_.get(j)));
//		}
//	}

//	public ArrayList<Pair<IdentificationVO,IdentificationVO>> getElutionOrderAssignments()
//	{
//		return elutionOrderAssignments_;
//	}
}
