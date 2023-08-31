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
import java.awt.Rectangle;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.analysis.exception.CalculationNotPossibleException;
import at.tugraz.genome.lda.utils.DoubleCalculator;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.vos.ResultAreaVO;
import at.tugraz.genome.lda.vos.ResultCompVO;
import at.tugraz.genome.lda.vos.ResultDisplaySettingsVO;
import at.tugraz.genome.maspectras.graphics.MacOSResizer;

/**
 * 
 * TODO: so basically I want to turn the heatmap gen upside down: generate each row (plus the names), then now knowing min and max, we initiate the gradient, paint it and then paint the rest
 * Rows should be their own objects that can be scrambled/sorted independently.
 * Compute median values for mol. species (this has to be independent of sum comp unfortunately), but we can use sum comp areas.
 * Same for C=C assigned stuff.
 * 
 * 
 * @author Juergen Hartler
 *
 */
public class LipidomicsHeatMap
{
	public final static int NUMBER_OF_COLOR_GROUPS = 9;
  private Color borderColor_=new Color(0,0,0);
  
  protected int defaultNrOfGradientPixels_ = 10000;
  
  protected final static Color DEFAULT_COLOR0 = new Color(000,040,180);
  protected final static Color DEFAULT_COLOR1 = new Color(000,153,102);
  protected final static Color DEFAULT_COLOR2 = new Color(255,255,000);
  protected final static Color DEFAULT_COLOR3 = new Color(255,255,000);
  protected final static Color DEFAULT_COLOR4 = new Color(255,000,000);
  private final static String LABEL_DB_CONFLICT = "\u03C9-DB assignment conflict";
  private final static String LABEL_DB_ASSIGNED = "\u03C9-DB assigned";
  
  protected int heatRectWidth_ = 20;
  protected int heatRectHeight_ = 10;
  protected int borderCorrection_ = 1;
  
  protected int annotationSpace_ = 200;
  protected int sampleSpace_ = 80;
  protected int pictureIndent_ = 5;
  protected int gradientHeight_ = 10;
  protected int analyteNameSpace_ = 3;
  
  protected int analyteTextWidthMax_ = 0;
  protected int sampleNameStart_ = 0;
  protected int doubleBondTextWidthMax_ = 0;
  
  int partCorrection_ = 0;
	
	protected Vector<String> sampleNames_;
  protected ArrayList<String> analyteNames_;
  
  protected BufferedImage gradient_;
  protected transient BufferedImage expressionImage_;
	
  protected float lowerValue_ = 0.2f;
  protected float upperValue_ = 5f;
  protected SampleLookup lookup_;
  private ArrayList<HeatMapRow> heatMapRows_;
  private boolean isMarkDoublePeaks_;
  
  protected boolean ignorePlatformSettings_;

  public final static Color ATTENTION_COLOR_DOUBLE_PEAK = Color.YELLOW;
  public final static Color ATTENTION_COLOR_NOT_ALL_MODS = new Color(51,153,255);
  public final static Color ATTENTION_DOUBLE_AND_NOT_ALL_MODS = new Color(204,153,255);
  
  public LipidomicsHeatMap(Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfOneGroup, 
  		Vector<String> experimentNames, SampleLookup lookup, Vector<String> moleculeNames,
  		int maxIsotope, ResultDisplaySettingsVO settingVO, boolean isMarkDoublePeaks) throws CalculationNotPossibleException
  {
  	this.isMarkDoublePeaks_ = isMarkDoublePeaks;
  	this.sampleNames_ = experimentNames;
  	this.heatMapRows_ = computeRows(resultsOfOneGroup, moleculeNames, maxIsotope, settingVO);
    this.analyteNames_ = computeAnalyteNames();
    this.gradient_ = this.createThreeColGradientImage(this.defaultNrOfGradientPixels_);
    this.lookup_ = lookup;
    this.ignorePlatformSettings_ = false;
  }
  
  private ArrayList<String> computeAnalyteNames()
  {
  	ArrayList<String> analyteNames = new ArrayList<String>();
  	for (HeatMapRow row : heatMapRows_)
    {
  		analyteNames.add(row.getAnalyteName());
    }
  	return analyteNames;
  }
  
  /**
   * Computes a row for the heat map.
   * @param resultsOfOneGroup
   * @param experimentNames
   * @param moleculeNames
   * @param maxIsotope
   * @param settingVO
   * @return
   * @throws CalculationNotPossibleException
   */
  private ArrayList<HeatMapRow> computeRows(Hashtable<String,Hashtable<String,ResultCompVO>> resultsOfOneGroup,
      Vector<String> moleculeNames, int maxIsotope, ResultDisplaySettingsVO settingVO) throws CalculationNotPossibleException
  {
  	ArrayList<HeatMapRow> rows = new ArrayList<HeatMapRow>();
    for (int j=0; j!=moleculeNames.size();j++)
    {
    	String sumCompositionName = moleculeNames.get(j);
    	Hashtable<String,ResultCompVO> resultsOfOneSumComp = resultsOfOneGroup.get(sumCompositionName);
    	Hashtable<String,ArrayList<Double>> medianValues = computeMedianAreaAtMolSpeciesLevel(resultsOfOneSumComp, settingVO, maxIsotope);
    	ResultCompVO dummyCompVO = new ResultCompVO();
    	
    	HeatMapRow rowSumComp = computeHeatMapRow(ResultCompVO.SUM_COMPOSITION, medianValues, sumCompositionName, settingVO, resultsOfOneSumComp, maxIsotope, dummyCompVO);
    	rows.add(rowSumComp);
    	
    	for (String name : medianValues.keySet())
    	{
    		if (name.equals(ResultCompVO.SUM_COMPOSITION)) continue;
    		HeatMapRow row = computeHeatMapRow(name, medianValues, sumCompositionName, settingVO, resultsOfOneSumComp, maxIsotope, dummyCompVO);
    		row.setAttentionValues(rowSumComp.getAttentionValues());
    		rows.add(row);
    	}
    }
    return rows;
  }
  
  private HeatMapRow computeHeatMapRow(String name, Hashtable<String,ArrayList<Double>> medianValues, String sumCompositionName, ResultDisplaySettingsVO settingVO, 
  		Hashtable<String,ResultCompVO> resultsOfOneSumComp, int maxIsotope, ResultCompVO dummyCompVO)
  {
  	ArrayList<Double> median = medianValues.get(name);
		HeatMapRow row = new HeatMapRow(sumCompositionName, name, computePreferredUnit(median.get(0), settingVO));
		Hashtable<String,Color> attentionValues = new Hashtable<String,Color>();
		for (int i=0;i!=this.sampleNames_.size();i++)
    {
			if (row.getSumCompositionName().equals("28:0") && row.getRtGroup().equals("12.83"))
			{
				System.out.println(row.getAnalyteName());
			}
			String experiment = this.sampleNames_.get(i);
			double relativeValue = -1d;
			if (isReasonableIdentification(name, resultsOfOneSumComp, experiment))
      {
      	ResultCompVO compVO = resultsOfOneSumComp.get(experiment);
      	compVO.addRelativeMedianArea(name, median);
      	row.addCompVO(experiment, compVO);
      	if (name.equals(ResultCompVO.SUM_COMPOSITION) && this.isMarkDoublePeaks_)
      	{
      		computeAttentionValuesRow(attentionValues, experiment, compVO, maxIsotope, sumCompositionName);
      	}
      	else
      	{
      		row.addMolecularSpeciesContributionOfAllMods(experiment, compVO.getResultMolecule().getMolecularSpeciesContributionOfAllMods(name));
      	}
      	relativeValue = computeRelativeValue(compVO, maxIsotope, settingVO, name, row.getMolecularSpeciesContribution(experiment));
      	row.addRelativeValue(experiment, relativeValue);
      }
			else
      {
				dummyCompVO.addRelativeMedianArea(name, median);
      	row.addCompVO(this.sampleNames_.get(i), dummyCompVO);
      	row.addRelativeValue(this.sampleNames_.get(i), relativeValue);
      }
    }
		if (name.equals(ResultCompVO.SUM_COMPOSITION)) row.setAttentionValues(attentionValues);
		return row;
  }
  
  private boolean isReasonableIdentification(String name, Hashtable<String,ResultCompVO> resultsOfOneSumComp, String experiment)
  {
  	if (resultsOfOneSumComp.containsKey(experiment))
  	{
  		if (name.equals(ResultCompVO.SUM_COMPOSITION))
    	{
    		return true;
    	}
  		else if (resultsOfOneSumComp.get(experiment).getResultMolecule() == null)
    	{
    		return false;
    	}
  		else if (resultsOfOneSumComp.get(experiment).getResultMolecule().getMolecularSpeciesContributionOfAllMods(name).isInfinite() ||
    			resultsOfOneSumComp.get(experiment).getResultMolecule().getMolecularSpeciesContributionOfAllMods(name).isNaN())
    	{
    		return false;
    	}
  		else
  		{
  			return true;
  		}
  	}
  	return false;
  }
  
  private Hashtable<String,ArrayList<Double>> computeMedianAreaAtMolSpeciesLevel(Hashtable<String,ResultCompVO> result, 
  		ResultDisplaySettingsVO settingVO, int maxIsotope) throws CalculationNotPossibleException
  {
  	Hashtable<String,ArrayList<Double>> medianValues = new Hashtable<String,ArrayList<Double>>();
  	if (result == null) return medianValues;
  	for (int k=0;k!=computeMaxApplicableIso(maxIsotope, settingVO, result);k++)
    {
  		Hashtable<String,ArrayList<Double>> areaList = computeAreaListOfIsotope(result, settingVO, k);
  		for (String name : areaList.keySet())
  		{
  			ArrayList<Double> areas = areaList.get(name);
  			if (!medianValues.containsKey(name)) medianValues.put(name, new ArrayList<Double>());
  			medianValues.get(name).add(DoubleCalculator.median(areas));
  		}
    }
  	return medianValues;
  }
  
  private Hashtable<String,ArrayList<Double>> computeAreaListOfIsotope(Hashtable<String,ResultCompVO> result, 
  		ResultDisplaySettingsVO settingVO, int k) throws CalculationNotPossibleException
  {
  	Hashtable<String,ArrayList<Double>> values = new Hashtable<String,ArrayList<Double>>();
  	values.put(ResultCompVO.SUM_COMPOSITION, new ArrayList<Double>());
  	for (int i=0;i!=this.sampleNames_.size();i++)
    {
      if (!result.containsKey(this.sampleNames_.get(i))) continue;
      ResultCompVO compVO = result.get(this.sampleNames_.get(i));
      double value = compVO.getArea(k, settingVO);
      if (value>0)
      {
        values.get(ResultCompVO.SUM_COMPOSITION).add(value);
        ResultAreaVO areaVO = compVO.getResultMolecule();
  			Set<String> names = areaVO.getAllMolecularSpeciesNamesHumanReadable();
  			for (String name : names)
  			{
  				if (!values.containsKey(name)) values.put(name, new ArrayList<Double>());
  				values.get(name).add(value * areaVO.getMolecularSpeciesContributionOfAllMods(name));
  			}
      } 
    }
  	return values;
  }
  
  private Double computeRelativeValue(ResultCompVO compVO, int maxIsotope, ResultDisplaySettingsVO settingVO, String molecularSpecies, double molecularSpeciesContribution)
  {
  	double relativeValue = compVO.getRelativeValue(compVO.getAvailableIsotopeNr(maxIsotope), settingVO, molecularSpecies, molecularSpeciesContribution);
    if (Double.isInfinite(relativeValue)||Double.isNaN(relativeValue))
      relativeValue = -1d;
  	return relativeValue;
  }
  
  private void computeAttentionValuesRow(Hashtable<String,Color> attentionValues, String experimentName, ResultCompVO compVO, int maxIsotope, String moleculeName)
  {
  	if (!(compVO.getMoreThanOnePeak(compVO.getAvailableIsotopeNr(maxIsotope)) || !compVO.hasAllModsFound())) return;
    Color color = null;
    if (compVO.getMoreThanOnePeak(compVO.getAvailableIsotopeNr(maxIsotope)) && !compVO.hasAllModsFound())
      color = ATTENTION_DOUBLE_AND_NOT_ALL_MODS;
    else if (compVO.getMoreThanOnePeak(compVO.getAvailableIsotopeNr(maxIsotope)))
      color = ATTENTION_COLOR_DOUBLE_PEAK;
    else if (!compVO.hasAllModsFound())
      color = ATTENTION_COLOR_NOT_ALL_MODS;
    attentionValues.put(experimentName, color);
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
  
  
  public BufferedImage createImage() {
    return this.createImage(null);
  }
  
  public BufferedImage createImage(Graphics2D extGraphics) 
  {
    BufferedImage dummyImage = new BufferedImage(1000,1000,BufferedImage.TYPE_3BYTE_BGR);
    Font descriptionFont = new Font("Dialog",Font.PLAIN, 9);
    Graphics2D dummyGraphics = (Graphics2D)dummyImage.getGraphics();
    dummyGraphics.setFont(descriptionFont);
    FontMetrics descriptionFontMetrics = dummyGraphics.getFontMetrics();
    int textWidth = 0;
    for (String sampleName : sampleNames_){
      String  toDisplay = getSampleNameToDisplay(sampleName);
      int width = descriptionFontMetrics.stringWidth(toDisplay);
      if (width>textWidth){
        textWidth = width;
      }  
    }
    int annotationWidth = 0;
    for (HeatMapRow row : heatMapRows_)
    {
    	int width = descriptionFontMetrics.stringWidth(row.getAnalyteName());
    	if (width>annotationWidth){
        annotationWidth = width;
      } 
    }
//    int doubleBondTextWidthMax = 0;
//    if (LipidParameterSet.isOmegaInformationAvailable()) {
//      doubleBondTextWidthMax = descriptionFontMetrics.stringWidth(LABEL_DB_CONFLICT);
//    }
    try 
    {
      this.sampleNameStart_ = this.pictureIndent_+this.gradientHeight_+descriptionFontMetrics.getHeight();
      this.sampleSpace_ = textWidth+sampleNameStart_;
      this.annotationSpace_ = annotationWidth+2*this.pictureIndent_;
      //      this.annotationSpace_ = annotationWidth+doubleBondTextWidthMax+2*this.pictureIndent_;
      Graphics2D expressionGraphics;
      int xcorrection=0;
      int ycorrection=0;
      int xSpaceForMap = this.heatRectWidth_*this.sampleNames_.size()+borderCorrection_;
      int imageSizeX = xSpaceForMap+xcorrection+annotationSpace_+2*pictureIndent_;
      expressionImage_=new BufferedImage(imageSizeX,
                                         this.heatRectHeight_*(this.heatMapRows_.size())+borderCorrection_+ycorrection+sampleSpace_+2*pictureIndent_,
                                         BufferedImage.TYPE_3BYTE_BGR);
      expressionGraphics=(Graphics2D)expressionImage_.getGraphics();
      if (extGraphics!=null)
      	expressionGraphics=extGraphics;
      expressionGraphics.setColor(getBackgroundColor());
      expressionGraphics.fillRect(0, 0, expressionImage_.getWidth(), expressionImage_.getHeight());
      paintHeader(expressionGraphics, imageSizeX-annotationSpace_);
      paintSampleNames(expressionGraphics,this.pictureIndent_ /*imageSizeX-pictureIndent_-xSpaceForMap-1*/,sampleSpace_);
      paintExpressionImage(expressionGraphics, this.pictureIndent_, sampleSpace_);
      paintAnalyteNames(expressionGraphics,this.pictureIndent_+xSpaceForMap,sampleSpace_);
      //       if (LipidParameterSet.isOmegaInformationAvailable()) {
      //         this.paintLipidDoubleBondNames(expressionGraphics,this.pictureIndent_+xSpaceForMap+analyteTextWidthMax_+doubleBondNameSpace_,sampleSpace_);
      //       }
    } catch (OutOfMemoryError e) { }
    return expressionImage_;      
  }
  
  protected void paintSampleNames(Graphics2D expressionGraphics, int xoffset, int sampleSpace_){
    Font bigFont = new Font("SansSerif", Font.PLAIN, 200);
    expressionGraphics.setFont(bigFont);
    FontMetrics fm2 = expressionGraphics.getFontMetrics();

    Font descriptionFont = new Font("Dialog",Font.PLAIN, 9);
    expressionGraphics.setColor(getLegendColor());
    
    expressionGraphics.setFont(descriptionFont);
    FontMetrics descriptionFontMetrics = expressionGraphics.getFontMetrics();
    int textHeight = descriptionFontMetrics.getHeight();
    double rotationAngle = -Math.PI/2.0;
    if (!isOSMacAndJavaLookAndFeel())
      expressionGraphics.rotate(rotationAngle);
    for (int i=0;i!=this.sampleNames_.size();i++){
      String sampleName = sampleNames_.get(i);
      String toDisplay = this.getSampleNameToDisplay(sampleName);
      if (isOSMacAndJavaLookAndFeel()){       
        BufferedImage image = new BufferedImage(fm2.stringWidth(toDisplay), fm2.getHeight(),
            BufferedImage.TYPE_INT_ARGB);
        Graphics2D g22 = (Graphics2D) image.getGraphics();
        g22.setColor(getLegendColor());
        g22.setFont(bigFont);
        g22.drawString(toDisplay,0,fm2.getAscent());
        int sw2 = fm2.stringWidth(toDisplay);
        int swXExtension2 = (int)(((double)sw2)*Math.cos(-rotationAngle));
        int swYExtension2 = (int)(((double)sw2)*Math.sin(-rotationAngle));
        int shYExtension2 = (int)(((double)fm2.getHeight())*Math.cos(-rotationAngle));
        int shXExtension2 = (int)(((double)fm2.getHeight())*Math.sin(-rotationAngle));
      
        BufferedImage image2 = new BufferedImage(swXExtension2+shXExtension2, swYExtension2+shYExtension2,
            BufferedImage.TYPE_INT_ARGB);
        Graphics2D g23 = (Graphics2D) image2.getGraphics();
        g23.rotate(rotationAngle);
        int normalStringXCoordinate2 = 0;
        int normalStringYCoordinate2 = swYExtension2;
        double coordinateRadius2 = Math.sqrt(Math.pow(0, 2)+Math.pow(normalStringYCoordinate2, 2));
        double angle2 = Math.atan((double)normalStringYCoordinate2/(double)normalStringXCoordinate2);
        g23.drawImage(image,(int)(coordinateRadius2*Math.cos(angle2-rotationAngle)),(int)(coordinateRadius2*Math.sin(angle2-rotationAngle)),null);
        g23.rotate(-rotationAngle);   
        //expressionGraphics.drawImage(image2.getScaledInstance(descriptionFontMetrics.getHeight(), descriptionFontMetrics.stringWidth(toDisplay), 0),this.heatRectWidth_ * i +this.heatRectWidth_/2,sampleSpace_-2-descriptionFontMetrics.stringWidth(toDisplay),null);
        image2 = MacOSResizer.resizeTrick(image2, descriptionFontMetrics.getHeight(), descriptionFontMetrics.stringWidth(toDisplay));
        expressionGraphics.drawImage(image2,this.heatRectWidth_ * i +this.heatRectWidth_/2,sampleSpace_-2-descriptionFontMetrics.stringWidth(toDisplay),null);
      }else{
        expressionGraphics.drawString(toDisplay,-sampleSpace_+2,this.heatRectWidth_ * i +this.heatRectWidth_/2+textHeight/2
            + xoffset);
      }
    }
    if (!isOSMacAndJavaLookAndFeel())
      expressionGraphics.rotate(-rotationAngle);
  }
  
  protected void paintExpressionImage(Graphics2D expressionGraphics, int xoffset, int yoffset)
  {
  	for (int i = 0; i < this.heatMapRows_.size(); i++) 
  	{
  		this.heatMapRows_.get(i).paintExpressionRow(expressionGraphics, xoffset, yoffset, i);
  	}
  }
  
  protected void paintAnalyteNames(Graphics2D expressionGraphics, int xoffset, int yoffset)
  {
    Font descriptionFont = new Font("Dialog",Font.PLAIN, 9);
    expressionGraphics.setColor(getLegendColor());
    expressionGraphics.setFont(descriptionFont);
    FontMetrics descriptionFontMetrics = expressionGraphics.getFontMetrics();
    int textHeight = descriptionFontMetrics.getHeight();
    for (int i=0; i<this.heatMapRows_.size(); i++)
    {
    	String analyteName = this.heatMapRows_.get(i).getAnalyteName();
    	int textWidth = descriptionFontMetrics.stringWidth(analyteName);
      if (textWidth>analyteTextWidthMax_)
        analyteTextWidthMax_ = textWidth;
      expressionGraphics.drawString(analyteName,xoffset+analyteNameSpace_,(this.heatRectHeight_*i) + 2 + yoffset+textHeight/2);
    }
  }
  
  public void paintAttentionRectangle(Graphics2D g2, Color color, int row, int column, int xoffset, int yoffset)
  {
  	this.heatMapRows_.get(row).paintAttentionRectangle(g2, color, row, column, xoffset, yoffset);
  }
  
  /**
   * returns the color of a heat map cell defined by row and column number
   * @param cellRow the row number of the cell
   * @param cellColumn the column number of the cell
   * @return color of a heat map cell defined by row and column number
   */
  public Color getColorForCell(int cellRow, int cellColumn){
    return new Color(getColorForValue(this.heatMapRows_.get(cellRow).getValue(sampleNames_.get(cellColumn))));
  }
  
  protected int getColorForValue(double value){
    int rgb1 = 0;
    if (value>1){
      int coordinate = getCoordinateBigger1(value);
      rgb1 = this.gradient_.getRGB(coordinate, 0);
    }else if (value>0){
      int coordinate = getCoordinateSmaller1(value);
      rgb1 = this.gradient_.getRGB(coordinate, 0);
    }else{
      rgb1 = Color.gray.getRGB();
    }
    return rgb1;
  }
  
  public int getExpressionImageXStart(){
    return this.pictureIndent_;
  }

  public int getExpressionImageXEnd(){
    return this.pictureIndent_+this.heatRectWidth_*sampleNames_.size();
  }
  
  public int getExpressionImageYStart(){
    return sampleSpace_;
  }

  public int getExpressionImageYEnd(){
    return this.sampleSpace_+this.heatRectHeight_*this.analyteNames_.size();
  }
  
  /**
   * Only to be used when hovering over analyte names!
   * @param x
   * @param y
   * @return
   */
  public String getRowName(int x, int y){
    String analyteName = null;
    if (this.getRowNameStart()<=x&&x<this.getRowNameEnd()&&
        this.getExpressionImageYStart()<=y&&y<this.getExpressionImageYEnd()){
    	analyteName = analyteNames_.get(this.calculateYPosition(y));
    }
    return analyteName;
  }
  
  /**
   * Only to be used when hovering over analyte names!
   * @param x
   * @param y
   * @return
   */
  public int getRowNumber(int x, int y)
  {
  	int rowNumber = -1;
  	if (this.getRowNameStart()<=x&&x<this.getRowNameEnd()&&
        this.getExpressionImageYStart()<=y&&y<this.getExpressionImageYEnd()){
  		rowNumber = this.calculateYPosition(y);
    }
  	return rowNumber;
  }
  
//  public Vector<String> getCellName(int x, int y){
//    Vector<String> cellName = new Vector<String>();
//    String sampleName = null;
//    String analyteName = null;
//    int[] cellPostion = this.getCellPosition(x, y);
//    if (cellPostion[0]>=0&&cellPostion[1]>=0){
//      sampleName = sampleNames_.get(cellPostion[0]);
//      analyteName = analyteNames_.get(cellPostion[1]);
//    }
//    cellName.add(sampleName);
//    cellName.add(analyteName);
//    return cellName;
//  }
  
  public int[] getCellPosition(int x, int y){
    int[] cellPosition = new int[2];
    cellPosition[0] = -1;
    cellPosition[1] = -1;
    if (this.getExpressionImageXStart()<=x&&x<this.getExpressionImageXEnd()&&
        this.getExpressionImageYStart()<=y&&y<this.getExpressionImageYEnd()){
      cellPosition[0] =this.calculateXPosition(x);
      cellPosition[1] = this.calculateYPosition(y);
    }
    return cellPosition;
  }
  
  public Rectangle getRectangleForCell(int x, int y){
    Rectangle rect = null;
    if (this.getExpressionImageXStart()<=x&&x<this.getExpressionImageXEnd()&&
        this.getExpressionImageYStart()<=y&&y<this.getExpressionImageYEnd()){
      int cellColumn = this.calculateXPosition(x);
      int cellRow = this.calculateYPosition(y);
      rect = getRectangleForCellByRowAndColumn(cellRow, cellColumn);
    }
    return rect;
  }
  
  /**
   * returns a rectangle of a heat map cell defined by row and column number
   * @param cellRow the row number of the cell
   * @param cellColumn the column number of the cell
   * @return rectangle of a heat map cell defined by row and column number
   */
  public Rectangle getRectangleForCellByRowAndColumn(int cellRow, int cellColumn){
    return new Rectangle(this.heatRectWidth_ * cellColumn+getExpressionImageXStart(),
        this.heatRectHeight_ * cellRow+getExpressionImageYStart(),this.heatRectWidth_,
        this.heatRectHeight_);
  }
  
  public int getRowNameStart(){
    return this.getExpressionImageXEnd()+this.analyteNameSpace_;
  }
  
  public int getRowNameEnd(){
    return this.getRowNameStart()+this.analyteTextWidthMax_;
  }
  
  public Rectangle getRectangleForRowName(String rowName){
    int rowPosition = -1;
    int count = 0;
    for (String name : this.analyteNames_){
      if (rowName.equalsIgnoreCase(name)){
        rowPosition = count;
        break;
      }
      count++;
    }
    Rectangle rect = null;
    if (rowPosition>=0){
      rect = new Rectangle(this.getRowNameStart(),this.heatRectHeight_ * rowPosition+getExpressionImageYStart(),
          this.analyteTextWidthMax_+this.analyteNameSpace_, this.heatRectHeight_);
    }
    return rect;
  }
  
  private int calculateXPosition(int x){
    return (x-this.getExpressionImageXStart())/heatRectWidth_;
  }  
  
  private int calculateYPosition(int y){
    return (y-this.getExpressionImageYStart())/heatRectHeight_;
  }

  public int getColumnNameStart(){
    return this.sampleNameStart_-2;
  }
  
  public int getColumnNameEnd(){
    return this.sampleSpace_-2;
  }

  public String getColumnName(int x, int y){
    String sampleName = null;
    if (this.getExpressionImageXStart()<=x&&x<this.getExpressionImageXEnd()&&
        this.getColumnNameStart()<=y&&y<this.getColumnNameEnd()){
      sampleName = sampleNames_.get(this.calculateXPosition(x));
    }
    return sampleName;
  }
  
  public Rectangle getRectangleForColumnName(String columnName){
    int columnPosition = -1;
    int count = 0;
    for (String name : this.sampleNames_){
      if (columnName.equalsIgnoreCase(name)){
        columnPosition = count;
        break;
      }
      count++;
    }
    Rectangle rect = null;
    if (columnPosition>=0){
      rect = new Rectangle(this.heatRectWidth_ * columnPosition+getExpressionImageXStart(),this.getColumnNameStart(),
          this.heatRectWidth_, this.getColumnNameEnd()-this.getColumnNameStart());
    }
    return rect;
  }
  
  private int getCoordinateBigger1(double intValue){
    int pixelsPositive = this.defaultNrOfGradientPixels_/2;
    int coordinate = pixelsPositive+(int)((intValue-1)*(pixelsPositive/(upperValue_-1)));
    if (coordinate>=this.defaultNrOfGradientPixels_||coordinate<0)
      coordinate = this.defaultNrOfGradientPixels_-1;
    return coordinate;
  }
  
  private int getCoordinateSmaller1(double intValue){
    int pixelsNegative = this.defaultNrOfGradientPixels_/2;
    intValue = 1/intValue;
    int coordinate = pixelsNegative-(int)((intValue-1)*(pixelsNegative/(1/lowerValue_-1)));
    if (coordinate<1)
      coordinate=1;
    return coordinate;
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
  
	public String getAbsoluteFilePathOfExperiment(int column)
	{
		for (HeatMapRow row : heatMapRows_)
		{
			ResultCompVO compVO = row.getCompVO(this.sampleNames_.get(column));
			if (compVO != null)
			{
				return compVO.getAbsoluteFilePath();
			}
		}
		return null;
	}
	
	public Color getAttentionProbe(int row, String experiment)
	{
		return this.heatMapRows_.get(row).getAttentionValues().get(experiment);
	}
	
	public Color getAttentionProbe(String analyte, String experiment)
	{
		return this.heatMapRows_.get(analyteNames_.indexOf(analyte)).getAttentionValues().get(experiment);
	}
	
	public ResultCompVO getCompVO(int column, int row)
	{
		return this.heatMapRows_.get(row).getCompVO(this.sampleNames_.get(column));
	}
	
	/**
	 * 
	 * @param name		the name of the experiment
	 * @param row
	 * @return
	 */
	public Double getMolecularSpeciesContributionOfAllMods(String name, int row)
	{
		return this.heatMapRows_.get(row).getMolecularSpeciesContribution(name);
	}
	
	public String getMolecularSpeciesLevelName(int row)
	{
		return this.heatMapRows_.get(row).getMolecularSpeciesName();
	}
	
	public String getSumCompositionName(int row)
	{
		return this.heatMapRows_.get(row).getSumCompositionName();
	}
	
	public String getRetentionTime(int row)
	{
		return this.heatMapRows_.get(row).getRtGroup();
	}
	
	public boolean isMolecularSpeciesLevel(int column, int row)
	{
		return !this.heatMapRows_.get(row).getMolecularSpeciesName().equals(ResultCompVO.SUM_COMPOSITION);
	}
	
	public String extractPreferredUnitForExp()
	{
    List<Integer> unitValues = new ArrayList<Integer>();
    for (HeatMapRow row : heatMapRows_)
    {
    	addUnitValue(unitValues, row.getPreferredUnit());
    }
    return computePreferredUnit(unitValues);
	}
	
	public String extractPreferredUnitForExp(Hashtable<String,String> preferredUnits)
	{
    List<Integer> unitValues = new ArrayList<Integer>();
    for (String unit : preferredUnits.values()){
    	addUnitValue(unitValues, unit);
    }
    return computePreferredUnit(unitValues);
  }
	
	private void addUnitValue(List<Integer> unitValues, String unit)
	{
		switch (unit)
  	{
  		case "":
  			unitValues.add(0);
  			break;
  		case "%":
  			unitValues.add(1);
  			break;
  		case "\u2030":
  			unitValues.add(2);
  			break;	
  		case "m":
  			unitValues.add(3);
  			break;
  		case "\u03BC":
  			unitValues.add(4);
  			break;
  		case "n":
  			unitValues.add(5);
  			break;
  		case "p":
  			unitValues.add(6);
  			break;
  		case "f":
  			unitValues.add(7);
  			break;
  		case "a":
  			unitValues.add(8);
  			break;
  		default:
  			break;
  	}
	}
	
	private String computePreferredUnit(List<Integer> unitValues)
	{
		String preferredUnit = "";
		Collections.sort(unitValues);
    int unitPlaceHolder = 0;
    if (unitValues.size()%2 == 0){
      unitPlaceHolder = unitValues.get(((unitValues.size()) / 2));
    }else{
      unitPlaceHolder = unitValues.get(((unitValues.size() + 1) / 2) - 1);
    }
    switch (unitPlaceHolder)
    {
    	case 0:
    		preferredUnit = "";
    		break;
    	case 1:
    		preferredUnit = "%";
    		break;
    	case 2:
    		preferredUnit = "\u2030";
    		break;
    	case 3:
    		preferredUnit = "m";
    		break;
    	case 4:
    		preferredUnit = "\u03BC";
    		break;
    	case 5:
    		preferredUnit = "n";
    		break;
    	case 6:
    		preferredUnit = "p";
    		break;
    	case 7:
    		preferredUnit = "f";
    		break;
    	case 8:
    		preferredUnit = "a";
    		break;
    	default:
    		break;
    }
    return preferredUnit;
	}
	
	public String getPreferredUnit(int row)
	{
		return this.heatMapRows_.get(row).getPreferredUnit();
	}
	
	public ArrayList<String> getAnalyteNames()
	{
		return analyteNames_;
	}
	
	public int getAnalyteIndex(String analyte)
	{
		return analyteNames_.indexOf(analyte);
	}
	
	public String getAnalyteName(int row)
	{
		return analyteNames_.get(row);
	}
	
	public String getOriginalAnalyteName(int row)
	{
		return heatMapRows_.get(row).getOriginalAnalyteName();
	}
	
	public String getExperimentName(int column)
	{
		String name = null;
		if (column >= 0 && column < sampleNames_.size()) name = sampleNames_.get(column);
		return name;
	}
	
	private class HeatMapRow
	{
		private String sumCompositionName_;
		private String rtGroup_;
		private String molecularSpeciesName_;
		private String preferredUnit_;
		private Hashtable<String,ResultCompVO> compVOs_; //expname to compVO
		private Hashtable<String,Color> attentionValues_; //expname to attentionvalue
		private Hashtable<String,Double> relativeValues_; //expname to rel value
		private Hashtable<String,Double> molecularSpeciesContribution_;
		private final static String RT_DELIMITER = "_";
		
		/**
		 * 
		 * @param moleculeName
		 * @param sumCompositionName
		 */
		private HeatMapRow(String sumCompositionName, String molecularSpeciesName, String preferredUnit)
		{
			this.sumCompositionName_ = extractSumCompositionName(sumCompositionName);
			this.rtGroup_ = extractRetentionTime(sumCompositionName);
			this.molecularSpeciesName_ = molecularSpeciesName;
			this.preferredUnit_ = preferredUnit;
			this.compVOs_ = new Hashtable<String,ResultCompVO>();
			this.attentionValues_ = new Hashtable<String,Color>();
			this.relativeValues_ = new Hashtable<String,Double>();
			this.molecularSpeciesContribution_ = new Hashtable<String,Double>();
		}
		
		private void paintExpressionRow(Graphics2D expressionGraphics, int xoffset, int yoffset, int row)
	  {
			double value;
	    expressionGraphics.setColor(borderColor_);
	    
	    for (int i = 0; i < sampleNames_.size(); i++)
	    {
	    	value = this.relativeValues_.get(sampleNames_.get(i));
	    	Color color = new Color(getColorForValue(value));
	    	expressionGraphics.setColor(color);
	    	expressionGraphics.fillRect(heatRectWidth_ * i + 1
            + xoffset + partCorrection_, heatRectHeight_ * row + 1
            + yoffset, heatRectWidth_ - 1, heatRectHeight_ - 1);
	    	if (this.attentionValues_.containsKey(sampleNames_.get(i))){
          paintAttentionRectangle(expressionGraphics,attentionValues_.get(sampleNames_.get(i)),row,i,xoffset,yoffset);
        } 
	    }
	  }
		
		/**
	   * this method paints an unfilled rectangle around the heat map rectangle containing the expression color
	   * @param g2 Graphics2 object for painting
	   * @param color the color to be used for painting the unfilled rectangle
	   * @param row the row number in the expression image
	   * @param column the column number in the expression image
	   * @param xoffset the offset in x direction of the heat map within the image
	   * @param yoffset the offset in y direction of the heat map within the image
	   */
	  private void paintAttentionRectangle(Graphics2D g2, Color color, int row, int column, int xoffset, int yoffset)
	  {
	    g2.setColor(this.attentionValues_.get(sampleNames_.get(column)));
	    g2.drawRect(heatRectWidth_ * column + 1
	        + xoffset + partCorrection_, heatRectHeight_ * row + 1
	        + yoffset, heatRectWidth_ - 2, heatRectHeight_ - 2);
	  }
	  
	  private String getAnalyteName()
	  {
	  	return String.format("%s ... %s min ... %s", getSumCompositionName(), getRtGroup(), getMolecularSpeciesName());
	  }
	  
	  private String getOriginalAnalyteName()
	  {
	  	return getSumCompositionName()+"_"+getRtGroup();
	  }
		
		private String extractSumCompositionName(String name)
		{
			return name.substring(0, name.indexOf(RT_DELIMITER));
		}
		
		private String extractRetentionTime(String name)
		{
			return name.substring(name.indexOf(RT_DELIMITER)+1, name.length());
		}

		private String getSumCompositionName()
		{
			return sumCompositionName_;
		}
		
		private String getMolecularSpeciesName()
		{
			return molecularSpeciesName_;
		}

		private String getRtGroup()
		{
			return rtGroup_;
		}

		private String getPreferredUnit()
		{
			return preferredUnit_;
		}

		private Hashtable<String,ResultCompVO> getCompVOs()
		{
			return compVOs_;
		}
		
		private void addCompVO(String experimentName, ResultCompVO compVO)
		{
			this.compVOs_.put(experimentName, compVO);
		}
		
		private ResultCompVO getCompVO(String experimentName)
		{
			return this.compVOs_.get(experimentName);
		}
		
//		private ArrayList<ResultCompVO> getCompVOsOfRow()
//		{
//			return new ArrayList<ResultCompVO>(compVOs_.values());
//		}

		private Hashtable<String,Color> getAttentionValues()
		{
			return attentionValues_;
		}

		private void setAttentionValues(Hashtable<String,Color> attentionValues)
		{
			if (this.molecularSpeciesName_.equals(ResultCompVO.SUM_COMPOSITION))
			{
				this.attentionValues_ = attentionValues;
			}
			else
			{
				for (String experiment : attentionValues.keySet())
				{
					if (molecularSpeciesContribution_.get(experiment) != null)
					{
						this.attentionValues_.put(experiment, attentionValues.get(experiment));
					}
				}
			}
		}

		private Hashtable<String,Double> getRelativeValues()
		{
			return relativeValues_;
		}
		
		private Double getValue(String experiment)
		{
			return relativeValues_.get(experiment);
		}
		
		private void addRelativeValue(String experimentName, Double relativeValue)
		{
			this.relativeValues_.put(experimentName, relativeValue);
		}
		
		private void addMolecularSpeciesContributionOfAllMods(String experimentName, Double value)
		{
			this.molecularSpeciesContribution_.put(experimentName, value);
		}
		
		private Double getMolecularSpeciesContribution(String experimentName)
		{
			if (this.molecularSpeciesName_.equals(ResultCompVO.SUM_COMPOSITION))
			{
				return 1.0;
			}
			else if (this.molecularSpeciesContribution_.get(experimentName) == null)
			{
				return 0.0;
			}
			return this.molecularSpeciesContribution_.get(experimentName);
		}
		
	}

}
