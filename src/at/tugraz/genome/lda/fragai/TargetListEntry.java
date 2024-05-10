package at.tugraz.genome.lda.fragai;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Objects;

public class TargetListEntry
{
	private String lipidClass_;
	private String species_;
	private Hashtable<String,Integer> sumFormula_;
	private ArrayList<Adduct> adducts_;
	private Double retentionTime_;
	private Double tolerance_;
	
	
	public TargetListEntry(String lipidClass, String species,
			Hashtable<String,Integer> sumFormula, ArrayList<Adduct> adducts,
			Double retentionTime, Double tolerance)
	{
		this.lipidClass_ = lipidClass;
		this.species_ = species;
		this.sumFormula_ = sumFormula;
		this.adducts_ = adducts;
		this.retentionTime_ = retentionTime;
		this.tolerance_ = tolerance;
	}


	public String getLipidClass()
	{
		return lipidClass_;
	}

	public String getSpecies()
	{
		return species_;
	}


	public Hashtable<String,Integer> getSumFormula()
	{
		return sumFormula_;
	}


	public ArrayList<Adduct> getAdducts()
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


	@Override
	public int hashCode()
	{
		return Objects.hash(adducts_, lipidClass_, species_, retentionTime_, sumFormula_,
				tolerance_);
	}


	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		TargetListEntry other = (TargetListEntry) obj;
		return Objects.equals(adducts_, other.adducts_)
				&& Objects.equals(lipidClass_, other.lipidClass_)
				&& Objects.equals(species_, other.species_)
				&& Objects.equals(retentionTime_, other.retentionTime_)
				&& Objects.equals(sumFormula_, other.sumFormula_)
				&& Objects.equals(tolerance_, other.tolerance_);
	}
	
}
