package at.tugraz.genome.lda.vos;

import java.util.Hashtable;
import java.util.Objects;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.utils.StaticUtils;

public class AdductVO
{
	private final static String ADDUCT_SEPARATOR = ";";
	private final static String ADDUCT_ADD = "+";
	private final static String ADDUCT_REMOVE = "-";
	private String adductName_;
	private String formulaString_;
	private int charge_;
	private Hashtable<String,Integer> formula_;
	
	public AdductVO(String name, String formula, int charge) throws ChemicalFormulaException
	{
		this.adductName_ = name;
		this.formulaString_ = formula;
		this.formula_ = StaticUtils.categorizeFormula(formula, true);
		this.charge_ = charge;
	}
	
	public AdductVO(String composite) throws ChemicalFormulaException
	{
		String[] split = composite.split(ADDUCT_SEPARATOR);
		this.adductName_ = split[0];
		this.formulaString_ = split[1];
		this.formula_ = StaticUtils.categorizeFormula(split[1], true);
		this.charge_ = Integer.parseInt(split[2]);
	}

	public String getAdductName()
	{
		return adductName_;
	}
	
	public String getFormulaString()
	{
		return formulaString_;
	}

	public int getCharge()
	{
		return charge_;
	}
	
	public Hashtable<String,Integer> getFormula()
	{
		return this.formula_;
	}

	@Override
	public int hashCode()
	{
		return Objects.hash(adductName_, charge_, formula_);
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
		AdductVO other = (AdductVO) obj;
		return Objects.equals(adductName_, other.adductName_)
				&& charge_ == other.charge_
				&& Objects.equals(formula_, other.formula_);
	}
	
	
}
