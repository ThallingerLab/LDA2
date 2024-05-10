package at.tugraz.genome.lda.fragai;

import java.util.Hashtable;
import java.util.Objects;

import at.tugraz.genome.lda.utils.StaticUtils;

public class Adduct
{
	private final static String ADDUCT_SEPARATOR = ";";
	private final static String ADDUCT_ADD = "+";
	private final static String ADDUCT_REMOVE = "-";
	private String adductName_;
	private int charge_;
	private Hashtable<String,Integer> addModifier_;
	private Hashtable<String,Integer> removeModifier_;
	
	public Adduct(String composite)
	{
		String[] split = composite.split(ADDUCT_SEPARATOR);
		this.adductName_ = split[0];
		String[] adductFormula = splitFormula(split[1]);
		this.charge_ = Integer.parseInt(split[2]);
		
		for (int i=0; i<adductFormula.length;i++)
		{
			if (adductFormula[i].startsWith(ADDUCT_ADD))
			{
				try {addModifier_ = StaticUtils.categorizeFormula(adductFormula[i].replace(ADDUCT_ADD,""));} catch (Exception ex) {}
			}
			else if (adductFormula[i].startsWith(ADDUCT_REMOVE))
			{
				try {removeModifier_ = StaticUtils.categorizeFormula(adductFormula[i].replace(ADDUCT_REMOVE,""));} catch (Exception ex) {}
			}
		}
	}
	
	public static String[] splitFormula(String formula) {
    int plusIndex = formula.indexOf(ADDUCT_ADD);
    int minusIndex = formula.indexOf(ADDUCT_REMOVE);

    String part1, part2;
    if (plusIndex != -1 && minusIndex != -1)
    {
    	if (plusIndex > minusIndex)
      {
      	part1 = formula.substring(minusIndex, plusIndex);
      	part2 = formula.substring(plusIndex);
      }
      else
      {
      	part1 = formula.substring(plusIndex, minusIndex);
      	part2 = formula.substring(minusIndex);
      }
    }
    else
    {
    	part1 = formula;
      part2 = "";
    }

    return new String[]{part1, part2};
	}

	public String getAdductName()
	{
		return adductName_;
	}

	public int getCharge()
	{
		return charge_;
	}

	public Hashtable<String,Integer> getAddModifier()
	{
		return addModifier_ == null ? new Hashtable<String,Integer>() : addModifier_;
	}

	public Hashtable<String,Integer> getRemoveModifier()
	{
		return removeModifier_ == null ? new Hashtable<String,Integer>() : removeModifier_;
	}

	@Override
	public int hashCode()
	{
		return Objects.hash(addModifier_, adductName_, charge_, removeModifier_);
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
		Adduct other = (Adduct) obj;
		return Objects.equals(addModifier_, other.addModifier_)
				&& Objects.equals(adductName_, other.adductName_)
				&& charge_ == other.charge_
				&& Objects.equals(removeModifier_, other.removeModifier_);
	}
	
	
}
