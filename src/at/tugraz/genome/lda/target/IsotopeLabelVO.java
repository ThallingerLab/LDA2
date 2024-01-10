package at.tugraz.genome.lda.target;

import java.util.Hashtable;
import java.util.LinkedHashMap;

import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * Class containing information about an isotope label
 * 
 * TODO: we do not need the labelElements_ because we have the labelId_, which is the prefix. 
 * Alternatively we could make it work without prefix, but the LDA kinda likes to differenciate between them visually, so we would have to assign a custom label... not worth it.
 * 
 * @author Leonida M. Lamp
 *
 */
public class IsotopeLabelVO
{
	/** the label elements including its amount for each label indicator*/
  protected Hashtable<String,Integer> labelElements_;
  /** the omega double bond position this label stands for*/  
  protected int omegaPosition_;
  /** the prefix that is used to indicate an isotope label*/
  protected String labelId_;
  
  /**
   * 
   * @param labelElements
   * @param omegaPosition
   * @param labelId
   */
  public IsotopeLabelVO(
  		Hashtable<String,Integer> labelElements, int omegaPosition, String labelId)
  {
  	this.labelElements_ = labelElements;
  	this.omegaPosition_ = omegaPosition;
  	this.labelId_ = labelId;
  }
  

  /**
   * getter for prefix that is used to indicate an isotope label
   * @return prefix that is used to indicate an isotope label
   */
  public String getLabelId()
  {
    return labelId_;
  }

  
  /**
   * setter for prefix that is used to indicate an isotope label
   * @param labelId prefix that is used to indicate an isotope label
   */
  public void setLabelId(String labelId)
  {
    this.labelId_ = labelId;
  }

  
  /**
   * getter for label elements including its amount for each label indicator
   * @return label elements including its amount for each label indicator
   */
  public Hashtable<String,Integer> getLabelElements()
  {
    return labelElements_;
  }

  
  /**
   * setter for label elements including its amount for each label indicator
   * @param labelElements label elements including its amount for each label indicator
   */
  public void setLabelElements(Hashtable<String,Integer> labelElements)
  {
    this.labelElements_ = labelElements;
  }
  
  
  /**
   * getter for the omega position this label stands for
   * @return omega position
   */
  public int getOmegaPosition()
  {
    return omegaPosition_;
  }

  
  /**
   * setter for the omega position this label stands for
   * @param omegaPosition the omega position
   */
  public void setOmegaPosition(int omegaPosition)
  {
    this.omegaPosition_ = omegaPosition;
  }

  
  /**
   * compares whether two value objects of this class are the same
   * @param other the other value object to compare to
   * @return true if both contain the same label and the same label elements
   */
  public boolean isEqual(IsotopeLabelVO other) {
    if (!labelId_.contentEquals(other.labelId_))
      return false;
    if (labelElements_.size()!=other.labelElements_.size())
      return false;
    for (String element : this.labelElements_.keySet()) {
      if (!other.labelElements_.containsKey(element))
        return false;
      if (this.labelElements_.get(element).intValue()!=other.labelElements_.get(element).intValue())
        return false;
    }
    return true;
  }
  
  
  public String toString() {
    return labelId_+": n-"+String.valueOf(omegaPosition_)+"; "+StaticUtils.getFormulaInHillNotation_PlusFirst(this.labelElements_, true);
  }
  
}
