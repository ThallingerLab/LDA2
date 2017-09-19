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

package at.tugraz.genome.lda.swing;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.math.BigDecimal;
import java.util.Collection;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JPanel;
import javax.swing.event.MouseInputAdapter;

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.analysis.SampleLookup;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class BarChart extends JPanel //implements ActionListener
{

  private static final long serialVersionUID = 5143957226504883251L;
  private boolean doubleSided_;
  private boolean logarithmic_;
  private double lowestValue_;
  private double highestValue_;
//  private Hashtable<String,Hashtable<String,Integer>> amountOfValues_;
  private Hashtable<String,Hashtable<String,Double>> valuesToPaint_;
  private Hashtable<String,Hashtable<String,Double>> sdsToPaint_;
  private String title_;
  private String subTitle_;
  private String normalizationText_;
  private Vector<String> molNames_;
  private Vector<String> originalValueNames_;
  SampleLookup sampleLookup_;
  private boolean subTitleHeightChanged_;
  
  private int leftMargin_ = 80;
  private int bottomMargin_ = 10;
  private int rightMargin_ = 15;
  private int topMargin_ = 30;
  private int xAxeDescrMargin_ = 5;
  private int roundValue_;
  
  private int marginLegend = 8;
  
  private String yAxisText_;
  private String sdText_;
  //private Hashtable<String,Color> colors_;
  ColorChooserDialog colors_;

  private int w0_;
  private int yLegendSize_;
  private int yLegendEntryHeight_;
  private int longestYLegend_;
  private int itemsYLegendRow_;
  
  private Rectangle rectToDraw_ = null;
  private int colorType_;
  protected boolean ignorePlatformSettings_;
  
  public BarChart(boolean doubleSided,boolean logarithmic,Hashtable<String,Hashtable<String,Double>> valuesToPaint, 
      Hashtable<String,Hashtable<String,Double>> sdsToPaint, String title, String subTitle, String normalizationText, 
      Vector<String> sortedMolNames, SampleLookup sampleLookup, Vector<String> sortedValueNames, double lowestTotalVale, double highestTotalValue, 
      //Hashtable<String,Hashtable<String,Integer>> amountOfValues, 
      String yAxisText, String sdText, ColorChooserDialog colorChooser_, int colorType){
    sampleLookup_ = sampleLookup;
    this.doubleSided_ = doubleSided;
    this.logarithmic_ = logarithmic;
    this.lowestValue_ = lowestTotalVale;
    this.highestValue_ = highestTotalValue;
//    amountOfValues_ = amountOfValues;
//    System.out.println("highest/lowest: "+highestValue_+";"+lowestValue_);
    title_ = title;
    subTitle_ = subTitle;
    normalizationText_ = normalizationText;
    molNames_ = sortedMolNames;
    valuesToPaint_ = valuesToPaint;
    sdsToPaint_ = sdsToPaint;
    colors_ = colorChooser_;
    sdText_ = sdText;
    if (sortedValueNames!=null){
      originalValueNames_ = sortedValueNames;
    }else{
      originalValueNames_ = new Vector<String>(valuesToPaint_.keySet());
      Collections.sort(originalValueNames_);
    }
    this.yAxisText_ = "[AU]";
    if (yAxisText!=null && yAxisText.length()>0)
      this.yAxisText_ = yAxisText;
    subTitleHeightChanged_ = false;
    MyListener myListener = new MyListener();
    this.addMouseListener(myListener);
    this.addMouseMotionListener(myListener);
    colorType_ = colorType;
    ignorePlatformSettings_ = false;
  }
  
  public void paint(Graphics g)
  {
    super.paint(g);
    this.drawBarchart(g);
  }
  
  public void drawBarchart(Graphics g, boolean ignorePlatformSpecificSettings){
    if  (ignorePlatformSpecificSettings) ignorePlatformSettings_ = true;
    drawBarchart(g);
    if  (ignorePlatformSpecificSettings) ignorePlatformSettings_ = false;
  }
  
  public void drawBarchart(Graphics g){
    Graphics2D g2 = (Graphics2D)g;
    g2.setColor(Color.WHITE);
    g2.fillRect(0, 0, this.getWidth(), this.getHeight());
    this.drawChart(g2);
  }
  
  public void drawChart(Graphics2D g2){
    Font bigFont = new Font("SansSerif", Font.PLAIN, 200);
    g2.setFont(bigFont);
    FontMetrics fm2 = g2.getFontMetrics();     
    g2.setFont(new Font("SansSerif", Font.PLAIN, 10));
    FontMetrics fmSmall = g2.getFontMetrics();
    
    g2.setColor(Color.BLACK);
    g2.setFont(new Font("SansSerif", Font.PLAIN, 12));
    FontMetrics fm = g2.getFontMetrics();
    int longestName = 0;
    Collection<String> xDescrNames = BarChart.extractDisplayNames(sampleLookup_,originalValueNames_);
    if (valuesToPaint_.size()>1) xDescrNames = molNames_;
    for (String name : xDescrNames){
      int sw = fm.stringWidth(name);
      if (sw>longestName)
        longestName = sw;
    }
    w0_ = this.getWidth()-leftMargin_-rightMargin_;
    int subTitleHeight = 0;
    if (subTitle_!=null && subTitle_.length()>0){
      g2.setFont(new Font("SansSerif", Font.PLAIN, 10));
      if (!subTitleHeightChanged_){     
        topMargin_ += fmSmall.getHeight()+1;
        subTitleHeightChanged_ = true;
      }
      subTitleHeight = fmSmall.getHeight()+1;
      g2.drawString(subTitle_,leftMargin_+w0_/2-fmSmall.stringWidth(subTitle_)/2,topMargin_-1-fmSmall.getHeight());
    }    
    g2.setFont(new Font("SansSerif", Font.BOLD, 14));
    FontMetrics fmBold = g2.getFontMetrics();
    g2.drawString(title_,leftMargin_+w0_/2-fmBold.stringWidth(title_)/2,topMargin_-1-fm.getHeight()-subTitleHeight);
    g2.setFont(new Font("SansSerif", Font.PLAIN, 10));
    int sizeForOnePart = w0_/originalValueNames_.size();
    int barDivisor = valuesToPaint_.size()*2+1;
    yLegendSize_ = 0;
    int coloredRectWidth = fm.getHeight()-4;
    int marginBetweenRectAndLegend = 4;
    yLegendEntryHeight_ = fm.getHeight();
    longestYLegend_ = 0;
    if (valuesToPaint_.size()>1){
      for (String valueName : BarChart.extractDisplayNames(sampleLookup_,originalValueNames_)){
        int sw = fm.stringWidth(valueName);
        if (sw>longestYLegend_)
          longestYLegend_ = sw; 
      }     
      longestYLegend_ +=coloredRectWidth+marginBetweenRectAndLegend+marginLegend;
      itemsYLegendRow_ = w0_/longestYLegend_;
      sizeForOnePart = w0_/valuesToPaint_.size();
      barDivisor = originalValueNames_.size()*2+1;
      // painting the legend
      int nrRows = originalValueNames_.size()/itemsYLegendRow_;
      if (originalValueNames_.size()%itemsYLegendRow_ != 0)
        nrRows++;
      yLegendSize_ = yLegendEntryHeight_*nrRows;
      for (int i=0; i!=originalValueNames_.size(); i++){
        int row = i/itemsYLegendRow_;
        int column = i%itemsYLegendRow_;
        Color myColor = colors_.getColor(colorType_, originalValueNames_.get(i));
        if (myColor.equals(Color.WHITE)){
          g2.drawRect(leftMargin_+longestYLegend_*column, this.getHeight()-bottomMargin_-(nrRows-row-1)*yLegendEntryHeight_-coloredRectWidth, coloredRectWidth, coloredRectWidth);
        }else{
          g2.setColor(colors_.getColor(colorType_, originalValueNames_.get(i)));
          g2.fillRect(leftMargin_+longestYLegend_*column, this.getHeight()-bottomMargin_-(nrRows-row-1)*yLegendEntryHeight_-coloredRectWidth, coloredRectWidth, coloredRectWidth);
          g2.setColor(Color.BLACK);
        }
        g2.drawString(BarChart.extractDisplayNames(sampleLookup_,originalValueNames_).get(i), leftMargin_+longestYLegend_*column+coloredRectWidth+marginBetweenRectAndLegend, this.getHeight()-bottomMargin_-(nrRows-1-row)*yLegendEntryHeight_);
      }
      
    }else{
      longestYLegend_ = fm.stringWidth(molNames_.get(0))+coloredRectWidth+marginBetweenRectAndLegend+marginLegend;
      itemsYLegendRow_ = 1;
      yLegendSize_ = yLegendEntryHeight_;
      if (colors_.getColor(ColorChooserDialog.DEFAULT_TYPE,"").equals(Color.WHITE))
        g2.drawRect(leftMargin_, this.getHeight()-bottomMargin_-coloredRectWidth, coloredRectWidth, coloredRectWidth);
      else{
        g2.setColor(colors_.getColor(ColorChooserDialog.DEFAULT_TYPE,""));
        g2.fillRect(leftMargin_, this.getHeight()-bottomMargin_-coloredRectWidth, coloredRectWidth, coloredRectWidth);
        g2.setColor(Color.BLACK);        
      }  
      g2.drawString(molNames_.get(0), leftMargin_+coloredRectWidth+marginBetweenRectAndLegend, this.getHeight()-bottomMargin_);
      
    }
    int sizeForOneBar = (sizeForOnePart*2)/barDivisor;
    double rotationAngle = -Math.PI/4;
//    if (doubleSided_)
//      rotationAngle = -Math.PI/6;
    int longestNameYExtension = (int)((double)longestName*Math.sin(-rotationAngle)+(double)fm.getHeight()*Math.cos(-rotationAngle));
    // horizontal axe
    int topMarginToUse = topMargin_;
    if (doubleSided_)
      topMarginToUse += (longestNameYExtension-5);    
    int h0 = this.getHeight()-topMarginToUse-bottomMargin_-longestNameYExtension-yLegendSize_-6;
    int y0 = 0;
    int rangeY = h0;
    
    if (normalizationText_!=null && normalizationText_.length()>0){
      double standRotAngle = -Math.PI/2;
      g2.setFont(new Font("SansSerif", Font.PLAIN, 10));
      int standX = fmSmall.getAscent()+5;
      int standY = topMarginToUse+h0/2+fmSmall.stringWidth(normalizationText_)/2;
      if (useMacSpecificSettings()){
        BufferedImage image = new BufferedImage(fm2.stringWidth(normalizationText_), fm2.getHeight(),
            BufferedImage.TYPE_INT_ARGB);
        Graphics2D g22 = (Graphics2D) image.getGraphics();
        g22.setColor(Color.BLACK);
        g22.setFont(bigFont);
        g22.drawString(normalizationText_,0,fm2.getAscent());
        int sw2 = fm2.stringWidth(normalizationText_);
        int swXExtension2 = (int)(((double)sw2)*Math.cos(-standRotAngle));
        int swYExtension2 = (int)(((double)sw2)*Math.sin(-standRotAngle));
        int shYExtension2 = (int)(((double)fm2.getHeight())*Math.cos(-standRotAngle));
        int shXExtension2 = (int)(((double)fm2.getHeight())*Math.sin(-standRotAngle));
      
        BufferedImage image2 = new BufferedImage(swXExtension2+shXExtension2, swYExtension2+shYExtension2,
            BufferedImage.TYPE_INT_ARGB);
        Graphics2D g23 = (Graphics2D) image2.getGraphics();
        g23.rotate(standRotAngle);
        int normalStringXCoordinate2 = 0;
        int normalStringYCoordinate2 = swYExtension2;
        double coordinateRadius2 = Math.sqrt(Math.pow(0, 2)+Math.pow(normalStringYCoordinate2, 2));
        double angle2 = Math.atan((double)normalStringYCoordinate2/(double)normalStringXCoordinate2);
        g23.drawImage(image,(int)(coordinateRadius2*Math.cos(angle2-standRotAngle)),(int)(coordinateRadius2*Math.sin(angle2-standRotAngle)),null);
        g23.rotate(-standRotAngle);
      
        g2.drawImage(image2.getScaledInstance(fmSmall.getHeight(), fmSmall.stringWidth(normalizationText_), 0),standX-fmSmall.getAscent(),standY-fmSmall.stringWidth(normalizationText_),null);
        
      }else{
        g2.rotate(standRotAngle);
        double standRadius = Math.sqrt(Math.pow(standX, 2)+Math.pow(standY, 2));
        double standAngle = Math.atan((double)standY/(double)standX);    
        g2.drawString(normalizationText_,(int)(standRadius*Math.cos(standAngle-standRotAngle)),(int)(standRadius*Math.sin(standAngle-standRotAngle)));
        g2.rotate(-standRotAngle);
      }
    }
    if (doubleSided_){
      y0 = h0/2+topMarginToUse;
      rangeY = h0/2;
    } else
      y0 = this.getHeight()-this.bottomMargin_-longestNameYExtension-yLegendSize_-6;
    g2.drawLine(this.leftMargin_-3, y0, this.leftMargin_ + w0_, y0);
    
    int yAxisYStart = this.getHeight()-bottomMargin_+3-longestNameYExtension-yLegendSize_-6;
    if (doubleSided_){
      yAxisYStart = yAxisYStart-3;
    }
    g2.drawLine(leftMargin_, yAxisYStart, leftMargin_, topMarginToUse);
    g2.drawLine(leftMargin_ - 2, topMarginToUse, leftMargin_ + 2, topMarginToUse);
    g2.drawLine(leftMargin_ - 2, topMarginToUse, leftMargin_, topMarginToUse-5);
    g2.drawLine(leftMargin_ + 2, topMarginToUse, leftMargin_, topMarginToUse-5);
    
    double in = 0;
    double din = 0;

    // **** Draw the y coordinates ****
    int x = leftMargin_;
    int j = fmSmall.stringWidth(yAxisText_);
    int yValue = y0 - h0 + 6;
    if (doubleSided_)
      yValue += h0/2;
    if (fmSmall.stringWidth(yAxisText_)>=x-4 || (sdText_!=null && sdText_.length()>0))
      yValue -= fmSmall.getHeight()-2;
    if (fmSmall.stringWidth(yAxisText_)>=x-4 && sdText_!=null && sdText_.length()>0 && doubleSided_)
      yValue -= fmSmall.getHeight();
    int axisDescrSpace = fmSmall.getHeight();
    if (fmSmall.stringWidth(yAxisText_)>=x-4){
      j = fmSmall.stringWidth(yAxisText_.substring(0,yAxisText_.lastIndexOf("[")-1));
      g2.drawString(yAxisText_.substring(0,yAxisText_.lastIndexOf("[")-1), x - j - 4, yValue);
      int unitWidth = fmSmall.stringWidth(yAxisText_.substring(yAxisText_.lastIndexOf("[")));
      g2.drawString(yAxisText_.substring(yAxisText_.lastIndexOf("[")), x - unitWidth - 4, yValue+fmSmall.getHeight());
      axisDescrSpace += fmSmall.getHeight();
    }else{
      g2.drawString(yAxisText_, x - j - 4, yValue);
    }
    if (sdText_!=null && sdText_.length()>0){
      j = fmSmall.stringWidth(sdText_);
      g2.drawString(sdText_, x - j - 4, yValue+axisDescrSpace);
    }
    int finalUsedYDescrSpace = yValue+axisDescrSpace+fmSmall.getHeight();
    
    double highestValueOfInterest = highestValue_;
    in = highestValueOfInterest / 10;
//    if (this.doubleSided_ && highestValue_<(1/this.lowestValue_)){
//      in = (1/this.lowestValue_)/5;
//    }
    
    if (doubleSided_){
      if (highestValue_<(1/this.lowestValue_))
        highestValueOfInterest = 1/this.lowestValue_;
      if (logarithmic_)
        in = Math.log10(highestValueOfInterest)/10;
      else
        in = (Math.log10(highestValueOfInterest)/Math.log10(2))/10;
    }
    roundValue_ = 2;
    if (in < 0.00000000000000000001){
      din = 0.00000000000000000001;
      roundValue_ = 20;
    }else if (in < 0.000000000000000000025){
      din = 0.000000000000000000025;
      roundValue_ = 21;
    } else if (in < 0.00000000000000000005){
      din = 0.00000000000000000005;
      roundValue_ = 20;
    } else if (in < 0.0000000000000000001){
      din = 0.0000000000000000001;
      roundValue_ = 19;
    }else if (in < 0.00000000000000000025){
      din = 0.00000000000000000025;
      roundValue_ = 20;
    } else if (in < 0.0000000000000000005){
      din = 0.0000000000000000005;
      roundValue_ = 19;
    } else if (in < 0.000000000000000001){
      din = 0.000000000000000001;
      roundValue_ = 18;
    }else if (in < 0.0000000000000000025){
      din = 0.0000000000000000025;
      roundValue_ = 19;
    } else if (in < 0.000000000000000005){
      din = 0.000000000000000005;
      roundValue_ = 18;
    } else if (in < 0.00000000000000001){
      din = 0.00000000000000001;
      roundValue_ = 17;
    }else if (in < 0.000000000000000025){
      din = 0.000000000000000025;
      roundValue_ = 18;
    } else if (in < 0.00000000000000005){
      din = 0.00000000000000005;
      roundValue_ = 17;
    } else if (in < 0.0000000000000001){
      din = 0.0000000000000001;
      roundValue_ = 16;
    }else if (in < 0.00000000000000025){
      din = 0.00000000000000025;
      roundValue_ = 17;
    } else if (in < 0.0000000000000005){
      din = 0.0000000000000005;
      roundValue_ = 16;
    } else if (in < 0.000000000000001){
      din = 0.000000000000001;
      roundValue_ = 15;
    }else if (in < 0.0000000000000025){
      din = 0.0000000000000025;
      roundValue_ = 16;
    } else if (in < 0.000000000000005){
      din = 0.000000000000005;
      roundValue_ = 15;
    } else if (in < 0.00000000000001){
      din = 0.00000000000001;
      roundValue_ = 14;
    }else if (in < 0.000000000000025){
      din = 0.000000000000025;
      roundValue_ = 15;
    } else if (in < 0.00000000000005){
      din = 0.00000000000005;
      roundValue_ = 14;
    } else if (in < 0.0000000000001){
      din = 0.0000000000001;
      roundValue_ = 13;
    }else if (in < 0.00000000000025){
      din = 0.00000000000025;
      roundValue_ = 14;
    } else if (in < 0.0000000000005){
      din = 0.0000000000005;
      roundValue_ = 13;
    } else if (in < 0.000000000001){
      din = 0.000000000001;
      roundValue_ = 12;
    }else if (in < 0.0000000000025){
      din = 0.0000000000025;
      roundValue_ = 13;
    } else if (in < 0.000000000005){
      din = 0.000000000005;
      roundValue_ = 12;
    } else if (in < 0.00000000001){
      din = 0.00000000001;
      roundValue_ = 11;
    }else if (in < 0.000000000025){
      din = 0.000000000025;
      roundValue_ = 12;
    } else if (in < 0.00000000005){
      din = 0.00000000005;
      roundValue_ = 11;
    } else if (in < 0.0000000001){
      din = 0.0000000001;
      roundValue_ = 10;
    }else if (in < 0.00000000025){
      din = 0.00000000025;
      roundValue_ = 11;
    } else if (in < 0.0000000005){
      din = 0.0000000005;
      roundValue_ = 10;
    } else if (in < 0.000000001){
      din = 0.000000001;
      roundValue_ = 9;
    }else if (in < 0.0000000025){
      din = 0.0000000025;
      roundValue_ = 10;
    } else if (in < 0.000000005){
      din = 0.000000005;
      roundValue_ = 9;
    } else if (in < 0.00000001){
      din = 0.00000001;
      roundValue_ = 8;
    }else if (in < 0.000000025){
      din = 0.000000025;
      roundValue_ = 9;
    } else if (in < 0.00000005){
      din = 0.00000005;
      roundValue_ = 8;
    }else if (in < 0.0000001){
      din = 0.0000001;
      roundValue_ = 7;
    }else if (in < 0.00000025){
      din = 0.00000025;
      roundValue_ = 8;
    } else if (in < 0.0000005){
      din = 0.0000005;
      roundValue_ = 7;
    } else if (in < 0.000001){
      din = 0.000001;
      roundValue_ = 6;
    }else if (in < 0.0000025){
      din = 0.0000025;
      roundValue_ = 7;
    } else if (in < 0.000005){
      din = 0.000005;
      roundValue_ = 6;
    } else if (in < 0.00001){
      din = 0.00001;
      roundValue_ = 5;
    }else if (in < 0.000025){
      din = 0.000025;
      roundValue_ = 6;
    } else if (in < 0.00005){
      din = 0.00005;
      roundValue_ = 5;
    } else if (in < 0.0001){
      din = 0.0001;
      roundValue_ = 4;
    }else if (in < 0.00025){
      din = 0.00025;
      roundValue_ = 5;
    } else if (in < 0.0005){
      din = 0.0005;
      roundValue_ = 4;
    } else if (in < 0.001){
      din = 0.001;
      roundValue_ = 3;
    }else if (in < 0.0025){
      din = 0.0025;
      roundValue_ = 4;
    } else if (in < 0.005){
      din = 0.005;
      roundValue_ = 3;
    }else if (in < 0.01)
      din = 0.01;
    else if (in < 0.025){
      din = 0.025;
      roundValue_ = 3;
    } else if (in < 0.05)
      din = 0.05;
    else if (in < 0.1)
      din = 0.1;
    else if (in < 0.25)
      din = 0.25;
    else if (in < 0.5)
      din = 0.5;
    else if (in < 1)
      din = 1;
    else if (in < 2.5)
      din = 2.5;
    else if (in < 5)
      din = 5;
    else if (in < 10)
      din = 10;
    else if (in < 25)
      din = 25;
    else if (in < 50)
      din = 50;
    else if (in < 100)
      din = 100;
    else if (in < 250)
      din = 250;
    else if (in < 500)
      din = 500;
    else if (in < 1000)
      din = 1000;
    else if (in < 2500)
      din = 2500;
    else if (in < 5000)
      din = 5000;
    else if (in < 10000)
      din = 10000;
    else if (in < 25000)
      din = 25000;
    else if (in < 50000)
      din = 50000;
    else if (in < 100000)
      din = 100000;
    else if (in < 250000)
      din = 250000;
    else if (in < 500000)
      din = 500000;
    else if (in < 1000000)
      din = 1000000;
    else if (in < 2500000)
      din = 2500000;
    else if (in < 5000000)
      din = 5000000;
    else if (in < 10000000)
      din = 10000000;
    else if (in < 25000000)
      din = 25000000;
    else if (in < 50000000)
      din = 50000000;
    else if (in < 100000000)
      din = 100000000;
    else if (in < 250000000)
      din = 250000000;
    else if (in < 500000000)
      din = 500000000;
    else if (in < 1000000000)
      din = 1000000000;
    else if (in < 2500000000l)
      din = 2500000000l;
    else if (in < 5000000000l)
      din = 5000000000l;
    else if (in < 10000000000l)
      din = 10000000000l;
    else if (in < 25000000000l)
      din = 25000000000l;
    else if (in < 50000000000l)
      din = 50000000000l;
    else if (in < 100000000000l)
      din = 100000000000l;
    else if (in < 250000000000l)
      din = 250000000000l;
    else if (in < 500000000000l)
      din = 500000000000l;
    else if (in < 1000000000000l)
      din = 1000000000000l;
    else if (in < 2500000000000l)
      din = 2500000000000l;
    else if (in < 5000000000000l)
      din = 5000000000000l;    
    else if (in < 10000000000000l)
      din = 10000000000000l;
    else if (in < 25000000000000l)
      din = 25000000000000l;
    else if (in < 50000000000000l)
      din = 50000000000000l;
    else if (in < 10000000000000l)
      din = 10000000000000l;
    else if (in < 25000000000000l)
      din = 25000000000000l;
    else if (in < 50000000000000l)
      din = 50000000000000l;    
    else if (in < 10000000000000l)
      din = 10000000000000l;
    else if (in < 25000000000000l)
      din = 25000000000000l;
    else if (in < 50000000000000l)
      din = 50000000000000l;
    else if (in < 10000000000000l)
      din = 10000000000000l;
    else if (in < 25000000000000l)
      din = 25000000000000l;
    else if (in < 50000000000000l)
      din = 50000000000000l;
    else if (in < 10000000000000l)
      din = 10000000000000l;
    else if (in < 25000000000000l)
      din = 25000000000000l;
    else if (in < 50000000000000l)
      din = 50000000000000l;
    else 
      din = 1000000000000000000l;
    if (logarithmic_&&!doubleSided_){
      din = 1;
    }

//    System.out.println("in/din: "+in+";"+din);
    int lowestExponent = getLowerOrHigherExponent(this.lowestValue_, 0, false);
    double valueRange = highestValueOfInterest-0;
    if (doubleSided_){
      valueRange = 0;
      if (logarithmic_){
        while (valueRange<Math.log10(highestValueOfInterest)){
          valueRange+=din;
        }
      }else{
        while (valueRange<(Math.log10(highestValueOfInterest)/Math.log10(2))){
          valueRange+=din;
        }
      }
//      if (logarithmic_)
//        highestExponent = getLowerOrHigherExponent(highestValueOfInterest, 0, true);
//      else      
//        highestExponent = getLowerOrHigherExponent(highestValueOfInterest, 0,2, true);
//      valueRange = highestExponent;
    }else if (logarithmic_){
      int highestExponent = getLowerOrHigherExponent(highestValueOfInterest, 0, true);
////      if (this.lowestValue_<1){
        ////valueRange = Math.pow(10,highestExponent-lowestExponent);
        valueRange = highestExponent-lowestExponent;       
////      }else{
////        valueRange = highestExponent;
////      }
//      System.out.println(this.highestValue_+";"+valueRange);
    }
    
    //else if (doubleSided_){
     // valueRange = M
   // }
    
    //else if (this.doubleSided_ && highestValue_<(1/this.lowestValue_)){
     // valueRange = 1/this.lowestValue_;
    //}
//    System.out.println("Value range: "+valueRange );
//    System.out.println(this.lowestValue_+";"+lowestExponent);
//    System.out.println("highestValueOfInterest: "+highestValueOfInterest);
//    System.out.println((new Double(Math.log10(lowestValue_))).);

    in = 0;
    int y = y0 - (int) (rangeY * in / valueRange);
//    if (doubleSided_){
//      in = 1;
//      y = y0 - (int) (rangeY * (in-1d) / (valueRange-1d));
//    }  
    if (logarithmic_)
      y = y0 - (int) (rangeY * in / valueRange);
    // while(in<m_maxIntensity)
    
    while (y > (y0 - rangeY)) {
      
      y = y0 - (int) ((rangeY * in) / valueRange);
      
//      if (logarithmic_)
//        y = y0 - (int) (rangeY * in /Math.log10(valueRange));
//      else if (doubleSided_)
//        y = y0 - (int) (rangeY * (in-1d) /(valueRange-1d));
      if (y < (y0 - rangeY + 10))
        break;
      if (in > Double.MAX_VALUE)
        break;

      g2.drawLine(x, y, x - 3, y);
      String s = "";
      if (doubleSided_){
        if (logarithmic_){
          s = "1.0E"+String.valueOf(Double.toString(roundDBL(in,roundValue_)));
        }else{
          s = "2^"+String.valueOf(Double.toString(roundDBL(in,roundValue_)));
        }
      }else if (logarithmic_ /*&& this.lowestValue_<1*/){
        s = "1.0E"+String.valueOf((int)in+lowestExponent);
      }
      else if (logarithmic_ && !doubleSided_)
        s = "1.0E"+String.valueOf((int)in);
      else if (in<10000000&&din>3){
        long inLong = (long)in;
        if (inLong<10000)
          s = Long.toString((long) in);
        else if (inLong<100000)
          s = Double.toString(roundDBL(in/10000,roundValue_))+"E4";
        else if (inLong<1000000)
          s = Double.toString(roundDBL(in/100000,roundValue_))+"E5";
        else
          s = Double.toString(roundDBL(in/1000000,roundValue_))+"E6";
//      else if (this.doubleSided_  && in<1.00000000001)
//        s = "1.0";
//      else if (this.doubleSided_ /*&& din< 2.5*/)
//        s = Double.toString(in);
      }else
        s = Double.toString(roundDBL(in,roundValue_));
//      System.out.println(in+";"+s+";"+(int) (rangeY * in / valueRange));
      int sw = fm.stringWidth(s);
      if (in == 0 ){
        if (y>finalUsedYDescrSpace)
          g2.drawString(s, x - sw - 4, y);
      } else {
        if ((y+4)>finalUsedYDescrSpace)
          g2.drawString(s, x - sw - 4, y + 4);
      }  
      in += din;
    }
    if (this.doubleSided_){
      in = -din;
      while (y < (y0 + rangeY)) {
        y = y0 - (int) ((rangeY * in) / valueRange);
//        if (logarithmic_)
//          y = y0 - (int) (rangeY * in /Math.log10(valueRange));
        if (y > (y0 + rangeY - 10))
          break;
        if (in > Double.MAX_VALUE)
          break;

        g2.drawLine(x, y, x - 3, y);
//        double inReverse = 1/(in);
//        if (din< 2.5)
//          inReverse = 1/(in+1);
////        int exponent = this.getLowerOrHigherExponent(inReverse, 0, false);
////        exponent--;
////        String doubleString = (new Double(Calculator.roundFloat((float)inReverse, -exponent))).toString();
//        int endIndex = doubleString.indexOf(".")-exponent+1;
//        if (endIndex>doubleString.length())
//          endIndex = doubleString.length();
//        String s = doubleString.substring(0,endIndex);
        String s = "";
        if (logarithmic_)
          s = "1.0E"+String.valueOf(Double.toString(roundDBL(in,roundValue_)));
        else
          s = "2^"+String.valueOf(Double.toString(roundDBL(in,roundValue_)));
        int sw = fm.stringWidth(s);
        if (y>finalUsedYDescrSpace)
          g2.drawString(s, x - sw - 4, y + 4);
        in -= din;      
      }
    }
    
    // painting the bars of the chart;
    int count = 0;
//    System.out.println("highestValue_: "+highestValue_);
    int amountOfMolecules = valuesToPaint_.size();
    
    int halfSDSize = sizeForOneBar/4;
    if (halfSDSize<1)
      halfSDSize = 1;
    for (String molName : molNames_){
      int expCount = 0;
      for (String name : originalValueNames_){
        double value = 0;
        double lowerSdValue = Double.NEGATIVE_INFINITY;
        double upperSdValue = Double.NEGATIVE_INFINITY;
        int ySDStart = -1;
        int ySDStop = -1;
        int xBarStart = this.leftMargin_+sizeForOnePart*count+sizeForOnePart/(barDivisor*2)+expCount*sizeForOneBar;
        if (this.doubleSided_){
          value = -1;
        }  
        if (valuesToPaint_.get(molName).containsKey(name)){
          value = valuesToPaint_.get(molName).get(name);
        }
        if (sdsToPaint_!=null && sdsToPaint_.get(molName).containsKey(name)){
          double sdValue = sdsToPaint_.get(molName).get(name);
          lowerSdValue = value - sdValue;
          upperSdValue = value + sdValue;
          
        
      }
//      double sdToPaint = -1;
//      if (upperSdValue>-1)
//        sdToPaint = sdValue;
      
        if (this.doubleSided_){
          double rightValue = value;
          if (logarithmic_){
            rightValue = Math.log10(rightValue);
            if (upperSdValue>Double.NEGATIVE_INFINITY){
              lowerSdValue = Math.log10(lowerSdValue);
              upperSdValue = Math.log10(upperSdValue);
            }
          }else{
            rightValue = Math.log10(rightValue)/Math.log10(2);
            if (upperSdValue>Double.NEGATIVE_INFINITY){
              lowerSdValue = Math.log10(lowerSdValue)/Math.log10(2);
              upperSdValue = Math.log10(upperSdValue)/Math.log10(2);
            }
          }
          int height = (int)((rangeY*rightValue)/valueRange);
          if (rightValue>0){
            drawFilledBar(g2,name,xBarStart, y0-height, sizeForOneBar, height);
          }else{
            drawFilledBar(g2,name,xBarStart, y0, sizeForOneBar, -height);
          }
          if (upperSdValue>Double.NEGATIVE_INFINITY){
            ySDStart = y0-(int)(((double)rangeY*(lowerSdValue))/valueRange);
//          if (ySDStart>y0)
//            ySDStart = y0;
            ySDStop = y0-(int)(((double)rangeY*(upperSdValue))/valueRange);

//          g2.drawLine(this.leftMargin_+sizeForOnePart*count+sizeForOnePart/2, ySDStart, this.leftMargin_+sizeForOnePart*count+sizeForOnePart/2, ySDStop);
//          g2.drawLine(this.leftMargin_+sizeForOnePart*count+sizeForOnePart/2-halfSDSize, ySDStart, this.leftMargin_+sizeForOnePart*count+sizeForOnePart/2+halfSDSize, ySDStart);
//          g2.drawLine(this.leftMargin_+sizeForOnePart*count+sizeForOnePart/2-halfSDSize, ySDStop, this.leftMargin_+sizeForOnePart*count+sizeForOnePart/2+halfSDSize, ySDStop);          
          }
//        System.out.println(lowerSdValue+";"+ySDStart);
//        System.out.println(upperSdValue+";"+ySDStop);
//        System.out.println(y0);

        }else{
          // this is for drawing the bar
          int height = (int)(((double)rangeY*value)/valueRange);
          if (logarithmic_){
          ////            if (lowestValue_<1){
              height = (int)((rangeY*Math.log10(value*Math.pow(10,-lowestExponent)))/valueRange); 
              if (upperSdValue>Double.NEGATIVE_INFINITY){
                lowerSdValue = Math.log10(lowerSdValue*Math.pow(10,-lowestExponent));
                upperSdValue = Math.log10(upperSdValue*Math.pow(10,-lowestExponent));
              ////              }  
////            }else{
            ////              height = (int)((rangeY*Math.log10(value))/valueRange);
            ////              if (upperSdValue>Double.NEGATIVE_INFINITY){
            ////                lowerSdValue = Math.log10(lowerSdValue);
            ////                upperSdValue = Math.log10(upperSdValue);
            ////              }
            }  
          }
          drawFilledBar(g2,name,xBarStart, y0-height, sizeForOneBar, height);
          // this is for get heights for SD
          if (upperSdValue>Double.NEGATIVE_INFINITY){
//          System.out.println("sdToPaint: "+sdToPaint);
            ySDStart = y0-(int)(((double)rangeY*(lowerSdValue))/valueRange);
            if (ySDStart>y0)
              ySDStart = y0;
            ySDStop = y0-(int)(((double)rangeY*(upperSdValue))/valueRange);
          }
        }
        if (upperSdValue>Double.NEGATIVE_INFINITY){
          g2.drawLine(xBarStart+sizeForOneBar/2, ySDStart, xBarStart+sizeForOneBar/2, ySDStop);
          g2.drawLine(xBarStart+sizeForOneBar/2-halfSDSize, ySDStart, xBarStart+sizeForOneBar/2+halfSDSize, ySDStart);
          g2.drawLine(xBarStart+sizeForOneBar/2-halfSDSize, ySDStop, xBarStart+sizeForOneBar/2+halfSDSize, ySDStop);
        }
        if (amountOfMolecules==1)
          count++;
        else
          expCount++;
      }
      if (amountOfMolecules>1)
        count++;
    }
    if (!useMacSpecificSettings())
      g2.rotate(rotationAngle);
    count = 0;
    
    for (String name : xDescrNames){
      int sw = fm.stringWidth(name);
      int swXExtension = (int)(((double)sw)*Math.cos(-rotationAngle));
      int swYExtension = (int)(((double)sw)*Math.sin(-rotationAngle));
      int shXExtension = (int)(((double)fm.getAscent())*Math.sin(-rotationAngle));
      int shYExtension = (int)(((double)fm.getAscent())*Math.cos(-rotationAngle));
      int lengthDifference = longestName-sw;
      int lengthDifferenceYExtension = (int)((double)(lengthDifference)*Math.sin(-rotationAngle));
      int normalStringXCoordinate = this.leftMargin_+sizeForOnePart*count+sizeForOnePart/2-swXExtension/2+shXExtension/2;
      int normalStringYCoordinate = this.getHeight()-bottomMargin_-lengthDifferenceYExtension/*+xAxeDescrMargin_*/-yLegendSize_;
      double value = 0;
      if (this.doubleSided_)
        value = -1;
      boolean lower1 = false;
//      boolean higher1 = false;
      double lowerValue = Double.MAX_VALUE;
      double upperValue = 0;
      if (valuesToPaint_.size()>1){
        for (String expName : originalValueNames_){
          if (valuesToPaint_.get(name).containsKey(expName)){
            double myValue = valuesToPaint_.get(name).get(expName);
            double upperMyValue = myValue;
            double lowerMyValue = myValue;
            if (myValue>0 && sdsToPaint_!=null && sdsToPaint_.get(name).containsKey(expName)){
              double sd = sdsToPaint_.get(name).get(expName);
              upperMyValue += sd;
              lowerMyValue -= sd;
            }  
            if (myValue>0 && lowerMyValue<lowerValue)
              lowerValue = lowerMyValue;
            if (upperMyValue > upperValue)
              upperValue = upperMyValue;
            if (this.doubleSided_){
              if (myValue>0 && myValue<1)
                lower1 = true; 
            }
          }
        }
        if (doubleSided_&&lower1&&upperValue>0&&(1/lowerValue>upperValue)){
          value = lowerValue;
        }else
          value = upperValue;
      }else{
        if (valuesToPaint_.values().iterator().next().containsKey(originalValueNames_.get(count)))
          value = valuesToPaint_.values().iterator().next().get(originalValueNames_.get(count));
      }
      
      // this is for the correction of the sample name position if sd is present
      if (sdsToPaint_!=null && sdsToPaint_.values().iterator().next().containsKey(name)){
        double sdValue = sdsToPaint_.values().iterator().next().get(name);
        if (value > 0 && !Double.isNaN(sdValue) && !Double.isInfinite(sdValue)){
          if (doubleSided_ && value<1)
            value -= sdValue;
          else
            value += sdValue;
        }
      }

      if (doubleSided_ && value > 0){
        double rightValue = value;
        if (logarithmic_)
          rightValue = Math.log10(rightValue);
        else{
          rightValue = Math.log10(rightValue)/Math.log10(2);
        }
        normalStringYCoordinate = y0-(int)((rangeY*rightValue)/valueRange)-xAxeDescrMargin_;
        if (value>0&&value<1)
          normalStringYCoordinate = y0-(int)((rangeY*rightValue)/valueRange)+2*xAxeDescrMargin_+swYExtension;
      }
      
      double coordinateRadius = Math.sqrt(Math.pow(normalStringXCoordinate, 2)+Math.pow(normalStringYCoordinate, 2));
      double angle = Math.atan((double)normalStringYCoordinate/(double)normalStringXCoordinate);
      if (useMacSpecificSettings()){
      
        BufferedImage image = new BufferedImage(fm2.stringWidth(name), fm2.getHeight(),
            BufferedImage.TYPE_INT_ARGB);
        Graphics2D g22 = (Graphics2D) image.getGraphics();
        g22.setColor(Color.BLACK);
        g22.setFont(bigFont);
        ////g22.setColor(Color.YELLOW);
        g22.drawString(name,0,fm2.getAscent());

        int sw2 = fm2.stringWidth(name);
        int swXExtension2 = (int)(((double)sw2)*Math.cos(-rotationAngle));
        int swYExtension2 = (int)(((double)sw2)*Math.sin(-rotationAngle));

        int shYExtension2 = (int)(((double)fm2.getHeight())*Math.cos(-rotationAngle));
        int shXExtension2 = (int)(((double)fm2.getHeight())*Math.sin(-rotationAngle));
      
        BufferedImage image2 = new BufferedImage(swXExtension2+shXExtension2, swYExtension2+shYExtension2,
            BufferedImage.TYPE_INT_ARGB);
        Graphics2D g23 = (Graphics2D) image2.getGraphics();
        ////g23.setColor(Color.RED);
        ////g23.fillRect(0, 0, swXExtension2+shXExtension2, swYExtension2+shYExtension2);
        g23.rotate(rotationAngle);
        int normalStringXCoordinate2 = 0;
        int normalStringYCoordinate2 = swYExtension2;
        double coordinateRadius2 = Math.sqrt(Math.pow(0, 2)+Math.pow(normalStringYCoordinate2, 2));
        double angle2 = Math.atan((double)normalStringYCoordinate2/(double)normalStringXCoordinate2);
        g23.drawImage(image,(int)(coordinateRadius2*Math.cos(angle2-rotationAngle)),(int)(coordinateRadius2*Math.sin(angle2-rotationAngle)),null);
        g23.rotate(-rotationAngle);
        g2.drawImage(image2.getScaledInstance(swXExtension+shXExtension, swYExtension+shYExtension, 0),normalStringXCoordinate-shXExtension,normalStringYCoordinate-swXExtension-shYExtension,null);
      }else
        g2.drawString(name,(int)(coordinateRadius*Math.cos(angle-rotationAngle)),(int)(coordinateRadius*Math.sin(angle-rotationAngle)));
      count++;
    }
    if (!useMacSpecificSettings())
      g2.rotate(-rotationAngle);
  }
  
  private void drawFilledBar(Graphics2D g2, String colName, int x, int y, int w, int h){
    Color color = Color.WHITE;
    if (valuesToPaint_.size()>1)
      color = colors_.getColor(colorType_, colName);
    else
      color = colors_.getColor(ColorChooserDialog.DEFAULT_TYPE,"");
    if (!color.equals(Color.WHITE)){
      g2.setColor(color);
      g2.fillRect(x, y, w, h);
      g2.setColor(Color.BLACK);      
    }
    g2.drawRect(x, y, w, h);
  }
  
  protected int getLowerOrHigherExponent(double value, int exponentDivFactor, boolean higher){
    return this.getLowerOrHigherExponent(value, exponentDivFactor,10, higher);
  }
  
  protected int getLowerOrHigherExponent(double value, int exponentDivFactor, double basis, boolean higher){
    int add = 0;
    if (higher)
      add++;
    if (value/basis>1) return this.getLowerOrHigherExponent(value/basis, ++exponentDivFactor,basis,higher);
    if (value/1>1){
      return exponentDivFactor+add;
    }else
      return this.getLowerOrHigherExponent(value*basis, --exponentDivFactor,basis,higher);
  }
  
  private double roundDBL(double targetDBL,int decimalPlace){
    BigDecimal bd = new BigDecimal(targetDBL);
    bd = bd.setScale(decimalPlace,BigDecimal.ROUND_HALF_UP);
    return (bd.doubleValue());
  }
  
  public BufferedImage getImage(){
    BufferedImage chart=new BufferedImage(getWidth(),getHeight(),
        BufferedImage.TYPE_3BYTE_BGR);
    drawBarchart(chart.getGraphics());
    return chart;
  }

  public int getRoundValue()
  {
    return roundValue_;
  }
  
  private class MyListener extends MouseInputAdapter {
    
    
    public void mouseMoved(MouseEvent e) {
      int x = e.getX();
      int y = e.getY();
      Graphics2D g2 = (Graphics2D)e.getComponent().getGraphics();
      if (rectToDraw_!=null){
        g2.setColor(Color.WHITE);
        g2.drawRect(rectToDraw_.x,rectToDraw_.y, rectToDraw_.width, rectToDraw_.height);
      }
      if (isInColorRegion(x,y)){
        int yInLegend = y-(getHeight()-bottomMargin_-yLegendSize_+2);
        int xInLegend = x-leftMargin_+marginLegend/2;
        int[] cellPosition = getYLegendCellPosition(xInLegend,yInLegend);
        String yLegendCellName = getCellNameForPosition(cellPosition);
        if (yLegendCellName!=null && yLegendCellName.length()>0){         
          rectToDraw_ = getRectangleForCell(cellPosition[0], cellPosition[1]);
          drawARectangle(g2);
        }
      }
    }
    
    public void mouseClicked(MouseEvent e) {
      int x = e.getX();
      int y = e.getY();
      if (isInColorRegion(x,y)){
        int yInLegend = y-(getHeight()-bottomMargin_-yLegendSize_+2);
        int xInLegend = x-leftMargin_+marginLegend/2;
        String yLegendCellName = getYLegendCellName(xInLegend,yInLegend);
        if (yLegendCellName!=null && yLegendCellName.length()>0){
          activateColorChooser(yLegendCellName);
        }
      }
    }   
  }
  
  private void activateColorChooser(String cellName){
    String title = "Choose your color for ";
    if (colorType_==ColorChooserDialog.DEFAULT_TYPE)
      title+="default";
    else if (colorType_==ColorChooserDialog.GROUP_TYPE)
      title+=cellName;
    else
      title+=sampleLookup_.getDisplayName(cellName);  
    colors_.showColorChooser(title, colorType_, cellName);
    //new ColorChooserDialog(new JFrame(),"Choose your color for "+cellName,initColor);
  }
  
  private void drawARectangle(Graphics2D g2){
    g2.setColor(Color.BLACK);
    g2.drawRect(rectToDraw_.x,rectToDraw_.y, rectToDraw_.width, rectToDraw_.height);      
  }
  
  private boolean isInColorRegion(int x, int y){
    boolean isThere = false;
    if (x>(leftMargin_-marginLegend/2) && x<(leftMargin_+w0_-marginLegend/2) &&
        y>(getHeight()-bottomMargin_-yLegendSize_+2)&&y<(getHeight()-bottomMargin_+2))
      isThere = true;
    return isThere;
  }
  
  private String getYLegendCellName(int x, int y){
    int[] cellPosition = getYLegendCellPosition(x,y);
    return getCellNameForPosition(cellPosition);
  }
  
  private String getCellNameForPosition(int[] cellPosition){
    String itemName = null;
    int amountOfItems = 1;
    if (valuesToPaint_.size()>1)
      amountOfItems = originalValueNames_.size();
    if (cellPosition[0]>-1 && cellPosition[0]<itemsYLegendRow_ && cellPosition[1]>-1 &&
        (cellPosition[1]*itemsYLegendRow_+cellPosition[0])<amountOfItems){
      if (valuesToPaint_.size()>1){
        itemName = originalValueNames_.get(cellPosition[1]*itemsYLegendRow_+cellPosition[0]);
      }else{
        itemName = molNames_.get(0);
      }
    }
    return itemName;    
  }
  
  /**
   * 
   * @param x
   * @param y
   * @return the position in the cell raster int[0]=column int[1]=row
   */
  private int[] getYLegendCellPosition(int x, int y){
    int[] positions = new int[2];
    positions[0] = getYLegendColumn(x);
    positions[1] = getYLegendRow(y);
    return positions;
  }
  
  private int getYLegendRow(int y){
    int row = y/yLegendEntryHeight_;
    return row;
  }
  
  private int getYLegendColumn(int x){
    int column = x/longestYLegend_;
    return column;
  }
  
  private Rectangle getRectangleForCell(int column, int row){
    Rectangle rect = new Rectangle(leftMargin_-marginLegend/2+column*longestYLegend_,
        getHeight()-bottomMargin_-yLegendSize_+2+yLegendEntryHeight_* row,longestYLegend_,
        yLegendEntryHeight_);
    return rect;
  }
  
  public static Vector<String> extractDisplayNames(SampleLookup sampleLookup, Vector<String> originalValueNames){
    Vector<String> displayNames = new Vector<String>();
    if (sampleLookup==null)
      displayNames = new Vector<String>(originalValueNames);
    else{
      for (String name: originalValueNames) displayNames.add(sampleLookup.getDisplayName(name));
    }
    return displayNames;
  }

//  public void actionPerformed(ActionEvent e)
//  {
//    if (e.getActionCommand().equalsIgnoreCase("AcceptColorSelection")){
//      colorDialog_.setVisible(false);
//      colorDialog_.dispose();
//    }else if(e.getActionCommand().equalsIgnoreCase("CancelColorSelection")){
//      colorDialog_.setVisible(false);
//      colorDialog_.dispose();
//    }
//  }
  
  private boolean useMacSpecificSettings(){
    if (!this.ignorePlatformSettings_&&Settings.isOSMacAndJavaLookAndFeel())
      return true;
    else return false;
  }
}
