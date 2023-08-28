/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome.lda.analysis;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.Hashtable;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ConcurrentHashMap;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.analysis.exception.CalculationNotPossibleException;
import at.tugraz.genome.lda.utils.DoubleCalculator;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;

/**
 * 
 * TODO: so basically I want to turn the heatmap gen upside down: generate each row (plus the names), then now knowing min and max, we initiate the gradient, paint it and then paint the rest
 * Rows should be their own objects that can be scrambled/sorted independently.
 * Compute median values for mol. species (this has to be independent of sum comp unfortunately), but we can use sum comp areas.
 * Same for C=C assigned stuff.
 * 
 * @author Juergen Hartler
 *
 */
public class LipidomicsHeatMap extends HeatMapGenerator
{
  protected float lowerValue_ = 0.2f;
  protected float upperValue_ = 5f;
  protected SampleLookup lookup_;
  private Hashtable<String,Hashtable<String,Color>> attentionProbes_;
  private ResultCompVO[][] compVOs_;
  private Hashtable<String,String> preferredUnit_;
  
  protected boolean ignorePlatformSettings_;

  public final static Color ATTENTION_COLOR_DOUBLE_PEAK = Color.YELLOW;
  public final static Color ATTENTION_COLOR_NOT_ALL_MODS = new Color(51,153,255);
  public final static Color ATTENTION_DOUBLE_AND_NOT_ALL_MODS = new Color(204,153,255);
  
  public LipidomicsHeatMap(Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfOneGroup, 
  		Vector<String> experimentNames, SampleLookup lookup, Vector<String> moleculeNames,
  		int maxIsotope, ResultDisplaySettingsVO settingVO, boolean isMarkDoublePeaks) throws CalculationNotPossibleException
  {
  	float[][] data = computeData(resultsOfOneGroup, experimentNames, moleculeNames, maxIsotope, settingVO);
  	setInputValues(data, experimentNames, moleculeNames, isMarkDoublePeaks ? attentionProbes_ : new Hashtable<String,Hashtable<String,Color>>());
    init();
    lookup_ = lookup;
    ignorePlatformSettings_ = false;
  }
  
  
  /**
   * Computes the data to be displayed and sets the values for attentionProbes_, compVOs_ and preferredUnit_.
   * @param resultsOfOneGroup
   * @param experimentNames
   * @param moleculeNames
   * @param maxIsotope
   * @param settingVO
   * @return
   * @throws CalculationNotPossibleException
   */
  private float[][] computeData(Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfOneGroup, Vector<String> experimentNames,
      Vector<String> moleculeNames, int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException
  {
  	float[][] data = new float[experimentNames.size()][moleculeNames.size()];
    Hashtable<String,Hashtable<String,Color>> attentionValues = new Hashtable<String,Hashtable<String,Color>>();
    Hashtable<String,String> preferredUnit = new Hashtable<String,String>();
    ResultCompVO[][] dataVOs = new ResultCompVO[experimentNames.size()][moleculeNames.size()];
    
    // this extracts the values for the median
    for (int j=0; j!=moleculeNames.size();j++)
    {
      Vector<Double> valuesForMedian = new Vector<Double>();
      if (moleculeNames.get(j).equals("32:2_14.41"))
      {
      	System.out.println("32:2_14.41");
      }
      Set<String> molecularSpeciesNames = findMolecularSpeciesNames(resultsOfOneGroup.get(moleculeNames.get(j)), experimentNames);
      //TODO: somewhere here I need to extract the mol. species names
      for (int k=0;k!=computeMaxApplicableIso(maxIsotope, settingVO, resultsOfOneGroup.get(moleculeNames.get(j)));k++)
      {
      	//TODO: here we gotta compute the median for molecular species as well
        double median = computeMedianOfResult(experimentNames, resultsOfOneGroup.get(moleculeNames.get(j)),settingVO, k);
        valuesForMedian.add(median);
        //TODO: this is fine if it applies to each sum. composition
        preferredUnit.put(moleculeNames.get(j), computePreferredUnit(median, settingVO));
      }
      
      for (int i=0;i!=experimentNames.size();i++)
      {
        double relativeValue = -1d;
        if (resultsOfOneGroup.containsKey(moleculeNames.get(j))&&
            resultsOfOneGroup.get(moleculeNames.get(j)).containsKey(experimentNames.get(i)))
        {
          ResultCompVO compVO = resultsOfOneGroup.get(moleculeNames.get(j)).get(experimentNames.get(i));
          if (compVO.getMoreThanOnePeak(compVO.getAvailableIsotopeNr(maxIsotope)) || !compVO.hasAllModsFound())
          {
          	computeAttentionValues(attentionValues, experimentNames.get(i), compVO, maxIsotope, moleculeNames.get(j));
          }
          //TODO: valuesForMedian needs to be added to compVO for each mol species
          compVO.setRelativeMedianArea(valuesForMedian);
          //TODO: the relativeValue needs to be computed for each mol. species as well
          relativeValue = computeRelativeValue(compVO, maxIsotope, settingVO);
          dataVOs[i][j] = compVO;
        }
        else
        {
          dataVOs[i][j] = new ResultCompVO();
          dataVOs[i][j].setRelativeMedianArea(valuesForMedian);
        }
        data[i][j] = (float)relativeValue;
      }
    }
    
    compVOs_ = dataVOs;
  	preferredUnit_ = preferredUnit;
  	attentionProbes_ = attentionValues;
    return data;
  }
  
  private Set<String> findMolecularSpeciesNames(Hashtable<String,ResultCompVO> result, Vector<String> experimentNames)
  {
  	Set<String> names = ConcurrentHashMap.newKeySet();
  	if (result != null)
  	{
  		for (int i=0;i!=experimentNames.size();i++)
      {
        if (result.containsKey(experimentNames.get(i)))
        {
          ResultCompVO compVO = result.get(experimentNames.get(i));
          ResultAreaVO areaVO = compVO.getResultMolecule();
          if (areaVO != null)
          {
          	Hashtable<String,Vector<Double>> rel = areaVO.computeMolecularSpeciesContributions();
          }
        }
      }
  	}
  	return names;
  }
  
  private Double computeRelativeValue(ResultCompVO compVO, int maxIsotope, ResultDisplaySettingsVO settingVO)
  {
  	double relativeValue = compVO.getRelativeValue(compVO.getAvailableIsotopeNr(maxIsotope), settingVO);
    if (Double.isInfinite(relativeValue)||Double.isNaN(relativeValue))
      relativeValue = -1d;
  	return relativeValue;
  }
  
  private void computeAttentionValues(Hashtable<String,Hashtable<String,Color>> attentionValues, String experimentName, ResultCompVO compVO, int maxIsotope, String moleculeName)
  {
  	Hashtable<String,Color> attention = new Hashtable<String,Color>();
    if (attentionValues.containsKey(experimentName))
      attention = attentionValues.get(experimentName);
    Color color = null;
    if (compVO.getMoreThanOnePeak(compVO.getAvailableIsotopeNr(maxIsotope)) && !compVO.hasAllModsFound())
      color = ATTENTION_DOUBLE_AND_NOT_ALL_MODS;
    else if (compVO.getMoreThanOnePeak(compVO.getAvailableIsotopeNr(maxIsotope)))
      color = ATTENTION_COLOR_DOUBLE_PEAK;
    else if (!compVO.hasAllModsFound())
      color = ATTENTION_COLOR_NOT_ALL_MODS;
    attention.put(moleculeName, color);
    attentionValues.put(experimentName, attention);
  }
  
  private Double computeMedianOfResult(Vector<String> experimentNames, Hashtable<String,ResultCompVO> result, 
  		ResultDisplaySettingsVO settingVO, int k) throws CalculationNotPossibleException
  {
  	Vector<Double> values = new Vector<Double>();
  	if (result != null)
  	{
  		for (int i=0;i!=experimentNames.size();i++)
      {
        if (result.containsKey(experimentNames.get(i)))
        {
        	//TODO: here, I could use a hashtable for the values - key sum.comp or mol.species name, values could then be added to each and returned
          ResultCompVO compVO = result.get(experimentNames.get(i));
          double standValue = compVO.getArea(k, settingVO);
          if (standValue>0)
          {
            values.add(standValue);
          }  
        }
      }
  	}
    return DoubleCalculator.median(values);
  }
  
  private int computeMaxApplicableIso(int maxIsotope, ResultDisplaySettingsVO settingVO, Hashtable<String,ResultCompVO> result)
  {
  	int maxApplicableIso = maxIsotope+1;
    if (!settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT) && 
        !settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT))
      maxApplicableIso = StaticUtils.getMaxApplicableIsotope(result, maxIsotope);
    return maxApplicableIso;
  }
  
  private String computePreferredUnit(double median, ResultDisplaySettingsVO settingVO)
  {
  	String unit = "";
    if (!(Double.isInfinite(median)||Double.isNaN(median)))
    {
      if (!(settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_VALUE) ||
          settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_BASE_PEAK)||
          settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT)||
          settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_HIGHEST_TOTAL_PEAK)||
          settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT)))
        unit = StaticUtils.extractPreferredUnit(median);
      else if (settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_BASE_PEAK)||
          settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_MEASURED_CLASS_AMOUNT)||
          settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_HIGHEST_TOTAL_PEAK)||
          settingVO.getType().equalsIgnoreCase(ResultDisplaySettingsVO.REL_TOTAL_AMOUNT))
        unit =  extractPercentUnit(median);
    }
    return unit;
  }
  
  
  protected BufferedImage createThreeColGradientImage(int sizeX) {
    
    BufferedImage gradient = new BufferedImage(sizeX, 1, BufferedImage.TYPE_3BYTE_BGR);
    Graphics2D graphics = gradient.createGraphics();
    GradientPaint gp0 = new GradientPaint(0, 0, Color.green, (sizeX)/2, 0,Color.BLACK);
    GradientPaint gp1 = new GradientPaint((sizeX)/2, 0, Color.black, sizeX, 0,Color.red);
    graphics.setPaint(gp0);
    graphics.drawRect(0, 0, (sizeX)/2, 1);
    graphics.setPaint(gp1);
    graphics.drawRect((sizeX)/2, 0, (sizeX)/2, 1);
    graphics.dispose();    
    return gradient;
  }
  
  protected void paintHeader(Graphics2D expressionGraphics,int imageSizeX){
    //paints the gradient
    GradientPaint gp0 = new GradientPaint(this.pictureIndent_, 0, Color.green, (imageSizeX)/2, 0,Color.black);
    GradientPaint gp1 = new GradientPaint((imageSizeX)/2, 0, Color.black,  imageSizeX-2*this.pictureIndent_, 0, Color.red);

    int gradientYPosition = this.pictureIndent_;
    expressionGraphics.setPaint(gp0);
    expressionGraphics.fillRect(this.pictureIndent_, gradientYPosition, (imageSizeX-2*this.pictureIndent_)/2, gradientHeight_);
    expressionGraphics.setPaint(gp1);
    expressionGraphics.fillRect(this.pictureIndent_+(imageSizeX-2*this.pictureIndent_)/2, gradientYPosition, (imageSizeX-2*this.pictureIndent_)/2, gradientHeight_);
    
    //paints the text of the gradient
    String minText = String.valueOf(this.lowerValue_);
    String middleText = "1.0";
    String maxText = String.valueOf(this.upperValue_);
    Font descriptionFont = new Font("Dialog",Font.PLAIN, 8);
    expressionGraphics.setColor(getLegendColor());
    expressionGraphics.setFont(descriptionFont);
    FontMetrics descriptionFontMetrics = expressionGraphics.getFontMetrics();
    int textHeight = descriptionFontMetrics.getHeight();
    int textWidth = descriptionFontMetrics.stringWidth(minText);
    int textYPosition = gradientYPosition+gradientHeight_+textHeight/2+2;
    expressionGraphics.drawString(minText,this.pictureIndent_,textYPosition);
    textWidth = descriptionFontMetrics.stringWidth(middleText);
    expressionGraphics.drawString(middleText,imageSizeX/2-textWidth/2,textYPosition);
    textWidth = descriptionFontMetrics.stringWidth(maxText);
    expressionGraphics.drawString(maxText,imageSizeX-this.pictureIndent_-textWidth,textYPosition);

//    String intermediate = "0.25";
//    textWidth = descriptionFontMetrics.stringWidth(intermediate);
//    int coordinate = this.getCoordinateSmaller1(Float.valueOf(intermediate));
//    expressionGraphics.drawString(intermediate,this.pictureIndent_+((imageSizeX-2*this.pictureIndent_)*coordinate)/this.defaultNrOfGradientPixels_-textWidth/2,textYPosition);
//    intermediate = "0.33";
//    textWidth = descriptionFontMetrics.stringWidth(intermediate);
//    coordinate = this.getCoordinateSmaller1(Float.valueOf(intermediate));
//    expressionGraphics.drawString(intermediate,this.pictureIndent_+((imageSizeX-2*this.pictureIndent_)*coordinate)/this.defaultNrOfGradientPixels_-textWidth/2,textYPosition);
//    intermediate = "0.5";
//    textWidth = descriptionFontMetrics.stringWidth(intermediate);
//    coordinate = this.getCoordinateSmaller1(Float.valueOf(intermediate));
//    expressionGraphics.drawString(intermediate,this.pictureIndent_+((imageSizeX-2*this.pictureIndent_)*coordinate)/this.defaultNrOfGradientPixels_-textWidth/2,textYPosition);
//    intermediate = "2.0";
//    textWidth = descriptionFontMetrics.stringWidth(intermediate);
//    coordinate = this.getCoordinateBigger1(Float.valueOf(intermediate));
//    expressionGraphics.drawString(intermediate,this.pictureIndent_+((imageSizeX-2*this.pictureIndent_)*coordinate)/this.defaultNrOfGradientPixels_-textWidth/2,textYPosition);
//    intermediate = "3.0";
//    textWidth = descriptionFontMetrics.stringWidth(intermediate);
//    coordinate = this.getCoordinateBigger1(Float.valueOf(intermediate));
//    expressionGraphics.drawString(intermediate,this.pictureIndent_+((imageSizeX-2*this.pictureIndent_)*coordinate)/this.defaultNrOfGradientPixels_-textWidth/2,textYPosition);
//    intermediate = "4.0";
//    textWidth = descriptionFontMetrics.stringWidth(intermediate);
//    coordinate = this.getCoordinateBigger1(Float.valueOf(intermediate));
//    expressionGraphics.drawString(intermediate,this.pictureIndent_+((imageSizeX-2*this.pictureIndent_)*coordinate)/this.defaultNrOfGradientPixels_-textWidth/2,textYPosition);    
  }
  
  protected int getColorForValue(float value){
    float intValue = value;
    int rgb1 = 0;
    if (value>1){
      int coordinate = getCoordinateBigger1(intValue);
      rgb1 = this.gradient_.getRGB(coordinate, 0);
    }else if (value>0){
      int coordinate = getCoordinateSmaller1(intValue);
      rgb1 = this.gradient_.getRGB(coordinate, 0);
    }else{
      rgb1 = Color.gray.getRGB();
    }
    return rgb1;
    //return this.mapColorToScale(rgb1);
  }
  
  private int getCoordinateBigger1(float intValue){
    int pixelsPositive = this.defaultNrOfGradientPixels_/2;
    int coordinate = pixelsPositive+(int)((intValue-1)*(pixelsPositive/(upperValue_-1)));
    if (coordinate>=this.defaultNrOfGradientPixels_||coordinate<0)
      coordinate = this.defaultNrOfGradientPixels_-1;
    return coordinate;
  }
  
  private int getCoordinateSmaller1(float intValue){
    int pixelsNegative = this.defaultNrOfGradientPixels_/2;
    intValue = 1/intValue;
    int coordinate = pixelsNegative-(int)((intValue-1)*(pixelsNegative/(1/lowerValue_-1)));
    if (coordinate<1)
      coordinate=1;
    return coordinate;
  }
  
  protected float getColorScale(){
    return (float)(10000-1)/(this.maxValue_-this.minValue_);
  }
  
  public String getSampleNameToDisplay(String sampleName){
    return this.lookup_.getDisplayName(sampleName);
  }
  
  private static String extractPercentUnit(double value){
    if (value<=0){
      return null;
    }else if (value*100>=1)
      return "%";
    else
      return "\u2030";
  }
  
  public BufferedImage createImage(Graphics2D extGraphics, boolean ignorePlatformSpecificSettings) {
    if  (ignorePlatformSpecificSettings) ignorePlatformSettings_ = true; //TODO: does this do anything at all!?
    BufferedImage image = createImage(extGraphics);
    if  (ignorePlatformSpecificSettings) ignorePlatformSettings_ = false;
    return image;
  }
  
  protected boolean isOSMacAndJavaLookAndFeel(){
    if (ignorePlatformSettings_) return false;
    return Settings.isOSMacAndJavaLookAndFeel();
  }
  
  protected Color getLegendColor(){
    return Color.WHITE;
  }
  
  protected Color getBackgroundColor(){
    return Color.BLACK;
  }
  
  public Hashtable<String,Hashtable<String,Color>> getAttentionProbes()
	{
		return attentionProbes_;
	}

	public ResultCompVO[][] getCompVOs()
	{
		return compVOs_;
	}

	public Hashtable<String,String> getPreferredUnit()
	{
		return preferredUnit_;
	}

}
