package at.tugraz.genome.lda.fragai;

import java.util.Hashtable;
import java.util.List;

public class TargetListEntry
{
	private String lipidClass_;
	private Hashtable<String,Integer> sumFormula_;
	private List<Adduct> adducts_;
	private Double retentionTime_;
	private Double tolerance_;
	
	
	public TargetListEntry(String lipidClass,
			Hashtable<String,Integer> sumFormula, List<Adduct> adducts,
			Double retentionTime, Double tolerance)
	{
		this.lipidClass_ = lipidClass;
		this.sumFormula_ = sumFormula;
		this.adducts_ = adducts;
		this.retentionTime_ = retentionTime;
		this.tolerance_ = tolerance;
	}


	public String getLipidClass()
	{
		return lipidClass_;
	}


	public Hashtable<String,Integer> getSumFormula()
	{
		return sumFormula_;
	}


	public List<Adduct> getAdducts()
	{
		return adducts_;
	}


	public Double getRetentionTime()
	{
		return retentionTime_;
	}


	public Double getTolerance()
	{
		return tolerance_;
	}
	
	
	
	
	
}
