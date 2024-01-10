package at.tugraz.genome.lda.target.calibration;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
import at.tugraz.genome.lda.quantification.LipidParameterSet;

public class IdentificationVO
{
	private ArrayList<File> files_ = new ArrayList<File>();
	private String lipidClass_;
	private String lipidSpecies_;
	private boolean isMSnAvailable_ = false;
	private ArrayList<Double> retentionTimes_ = new ArrayList<Double>();
	private ArrayList<Double> areas_ = new ArrayList<Double>();
  //TODO: this should be a setting
	private double groupingParameter_ = 0.05;
	
	/**
	 * 
	 * @param file					the file this analyte stems from for sanity check purposes
	 * @param lipidClass
	 * @param param
	 */
	protected IdentificationVO(File file, String lipidClass, LipidParameterSet param)
	{
		this.files_.add(file);
		this.lipidClass_ = lipidClass;
		this.lipidSpecies_ = param.getNameStringWithoutRt();
		this.isMSnAvailable_ = param instanceof LipidomicsMSnSet;
		this.retentionTimes_.add(param.getPreciseRT());
		this.areas_.add(new Double(param.getArea()));
	}
	
	/**
	 * Adds an analyte to the object if its retention time is within the grouping parameter
	 * Make sure the lipid class and lipid species match!
	 * @param file			the file this analyte stems from for sanity check purposes
	 * @param param
	 * @return
	 */
	protected boolean addParam(File file, LipidParameterSet param)
	{
		Double paramRT = param.getPreciseRT();
		for (int i=0; i<retentionTimes_.size(); i++)
		{
			Double rT = retentionTimes_.get(i);
			if (paramRT >= (rT - groupingParameter_) &&
					paramRT <= (rT + groupingParameter_))
			{
				this.files_.add(file);
				this.isMSnAvailable_ = param instanceof LipidomicsMSnSet ? true : isMSnAvailable_;
				this.retentionTimes_.add(param.getPreciseRT());
	  		this.areas_.add(new Double(param.getArea()));
	  		return true;
			}
		}
		return false;
	}
	
	Double getAverageRT()
	{
		Double sum = 0.0;
		for (Double rT : retentionTimes_)
		{
			sum += rT;
		}
		return sum/retentionTimes_.size();
	}
	
	Double getAverageArea()
	{
		Double sum = 0.0;
		for (Double area : areas_)
		{
			sum += area;
		}
		return sum/areas_.size();
	}
	
	boolean isMSnAvailable()
	{
		return isMSnAvailable_;
	}
	
	/**
	 * If an identification was contributed to more than once by the same result file, 
	 * either the grouping parameter or the experiment is not suitable.
	 * @return 		true if the identification was contributed to only once per file
	 */
	protected boolean isIdentificationUnambiguous()
	{
		Set<File> uniqueFiles = new HashSet<File>(files_);
		return uniqueFiles.size() == files_.size();
	}
	
	/**
	 * Returns the number of files this identification was found in.
	 * @return
	 */
	int getIdentificationCount()
	{
		Set<File> uniqueFiles = new HashSet<File>(files_);
		return uniqueFiles.size();
	}
	
	protected String getLipidClass()
	{
		return this.lipidClass_;
	}
	
	protected String getLipidSpecies()
	{
		return this.lipidSpecies_;
	}

	@Override
	public String toString()
	{
		return "IdentificationVO [files_=" + files_ + ", lipidClass_="
				+ lipidClass_ + ", lipidSpecies_=" + lipidSpecies_
				+ ", isMSnAvailable_=" + isMSnAvailable_ + ", retentionTimes_="
				+ retentionTimes_ + ", areas_=" + areas_ + "]";
	}
}
