package at.tugraz.genome.lda.target.calibration;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import at.tugraz.genome.lda.utils.Pair;

public class MatchedIdentificationVO
{
	ArrayList<IdentificationVO> originals_;
	int maxOriginals_;
	ArrayList<IdentificationVO> matches_;
	int maxMatches_;
	Pair<IdentificationVO,IdentificationVO> highestConfidencePair_;
	int confidence_ = 0;
	ArrayList<Pair<IdentificationVO,IdentificationVO>> elutionOrderAssignments_ = new ArrayList<Pair<IdentificationVO,IdentificationVO>>();
	
	protected MatchedIdentificationVO(ArrayList<IdentificationVO> originals, int maxOriginals, ArrayList<IdentificationVO> matches, int maxMatches)
	{
		this.originals_ = originals;
		this.maxOriginals_ = maxOriginals;
		this.matches_ = matches;
		this.maxMatches_ = maxMatches;
		findHighestConfidencePair();
		findElutionOrderAssignments();
	}
	
	/**
	 * Has to be executed after the highest confidence pair is found.
	 * Identifications are sorted based on retention time. Then they are added to an ArrayList
	 * first descending and then ascending starting from the highest confidence pair.
	 */
	private void findElutionOrderAssignments()
	{
		//TODO: make sure the sorting is done in the correct order!
		Collections.sort(originals_, obtainRTComparator());
		Collections.sort(matches_, obtainRTComparator());
		int indexOriginal = originals_.indexOf(highestConfidencePair_.getKey());
		int indexMatch = matches_.indexOf(highestConfidencePair_.getValue());
		
		int j = indexMatch;
		for (int i = indexOriginal-1; i>=0; i--)
		{
			j--;
			if (j < 0) continue;
			elutionOrderAssignments_.add(new Pair<IdentificationVO,IdentificationVO>(originals_.get(i),matches_.get(j)));
		}
		j = indexMatch;
		for (int i = indexOriginal+1; i<originals_.size(); i++)
		{
			j++;
			if (j >= matches_.size()) continue;
			elutionOrderAssignments_.add(new Pair<IdentificationVO,IdentificationVO>(originals_.get(i),matches_.get(j)));
		}
	}
	
	private void findHighestConfidencePair()
	{
		Collections.sort(this.originals_, obtainAreaComparator());
		Collections.sort(this.matches_, obtainAreaComparator());
		IdentificationVO highestIntensityOriginal = originals_.get(0);
		IdentificationVO highestIntensityMatch = matches_.get(0);
		this.highestConfidencePair_ = new Pair<IdentificationVO,IdentificationVO>(highestIntensityOriginal, highestIntensityMatch);
		
		if (isAlwaysFound(highestIntensityOriginal, maxOriginals_) && 
				isAlwaysFound(highestIntensityMatch, maxMatches_) &&
				isIntensityElutionOrderCorrelated())
		{
			this.confidence_ = 2;
		}
		else
		{
			this.confidence_ = 0; //TODO: delete after and lower other to 1
		}
//		else if (isAlwaysFound(highestIntensityOriginal, maxOriginals_) && 
//				isAlwaysFound(highestIntensityMatch, maxMatches_) &&
//				isIntensityElutionOrderCorrelated()
//		{
//			this.confidence_ = 1;
//		}
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
	private Comparator<IdentificationVO> obtainRTComparator()
	{
		return new Comparator<IdentificationVO>() 
		{
			@Override
			public int compare(IdentificationVO o1, IdentificationVO o2)
			{
				return o1.getAverageRT().compareTo(o2.getAverageRT());
			}
		};
	}

	public Pair<IdentificationVO,IdentificationVO> getHighestConfidencePair()
	{
		return highestConfidencePair_;
	}

	public int getConfidence()
	{
		return confidence_;
	}

	public ArrayList<Pair<IdentificationVO,IdentificationVO>> getElutionOrderAssignments()
	{
		return elutionOrderAssignments_;
	}
}
