/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Hashtable;
import java.util.Vector;

import javax.imageio.ImageIO;
import javax.swing.JFrame;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.parser.LDAResultReader;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.lda.swing.ExportPanel;
import at.tugraz.genome.lda.swing.LipidomicsTableCellRenderer;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.chromaviewer.MSMapViewer;
import at.tugraz.genome.maspectras.graphics.MacOSResizer;
import at.tugraz.genome.maspectras.quantification.Analyzer;
import at.tugraz.genome.maspectras.quantification.CgChromatogram;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.ChromatogramReader;
import at.tugraz.genome.maspectras.utils.StringUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ChromExportThread extends Thread
{
  
  private String lipidClass_;
  private File fileToStore_;
  private Dimension dim_; 
  private Vector<String> expsToExport_;
  private Vector<String> resultsToUse_;
  private Vector<String> chromsToUse_;
  private Vector<String> analsToExport_;
  private String picType_;
  private boolean analInColumn_;
  
  private String errorString_ = null;
  private boolean finished_ = false;
  private int currentLipidCount_;
  private String currentLipid_;
  private String currentExperiment_;
  private int totalAmountOfLipids_;
  private boolean interrupt_;
  private Double rtTolerance_;
  private Hashtable<String,Integer> isLookup_;
  private Hashtable<String,Integer> esLookup_;
  
  
  public ChromExportThread(String lipidClass,File fileToStore,Dimension dim, Vector<String> expsToExport,Vector<String> resultsToUse, Vector<String> chromsToUse, Vector<String> analsToExport, String picType, boolean analInColumn, Double rtTolerance,
      Hashtable<String,Integer> isLookup, Hashtable<String,Integer> esLookup){
    lipidClass_ = lipidClass;
    fileToStore_ = fileToStore;
    dim_ = dim; 
    expsToExport_ = expsToExport;
    resultsToUse_ = resultsToUse;
    chromsToUse_ = chromsToUse;
    analsToExport_ = analsToExport;
    picType_ = picType;
    analInColumn_ = analInColumn;
    errorString_ = null;
    finished_ = false;
    rtTolerance_ = rtTolerance;
    isLookup_ = isLookup;
    esLookup_ = esLookup;
  }
  
  public boolean finished(){
    return this.finished_;
  }
  
  public String getErrorString(){
    return this.errorString_;
  }
  
  public int getCurrentLipidCount()
  {
    return currentLipidCount_;
  }

  public String getCurrentLipid()
  {
    return currentLipid_;
  }

  public String getCurrentExperiment()
  {
    return currentExperiment_;
  }

  public int getTotalAmountOfLipids()
  {
    return totalAmountOfLipids_;
  }

  public void run(){
    try {
      exportChromatograms(fileToStore_,dim_, expsToExport_, resultsToUse_, chromsToUse_, analsToExport_, picType_, analInColumn_);
    }
    catch (IOException iox) {
      errorString_ = iox.getMessage();
    }
    this.finished_ = true;
  }
  
  private void exportChromatograms(File fileToStore,Dimension dim, Vector<String> expsToExport,Vector<String> resultsToUse, Vector<String> chromsToUse, Vector<String> analsToExport, String picType, boolean analInColumn)throws IOException{
    currentLipidCount_ = 0;
    currentLipid_ = "";
    currentExperiment_ = "";
    totalAmountOfLipids_ = expsToExport.size()*analsToExport.size();    
    
    if (picType.equalsIgnoreCase(ExportPanel.EXPORT_PNG)){
      BufferedImage chroms = new BufferedImage(dim.width,dim.height,BufferedImage.TYPE_3BYTE_BGR);
      drawChromatograms(chroms.getGraphics(),false,dim,expsToExport,resultsToUse,chromsToUse,analsToExport,analInColumn);
      ImageIO.write(chroms, "PNG", fileToStore);
    }else if (picType.equalsIgnoreCase(ExportPanel.EXPORT_SVG)){
      DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
      // Create an instance of org.w3c.dom.Document.
      String svgNS = "http://www.w3.org/2000/svg";
      Document document = domImpl.createDocument(svgNS, "svg", null);
      // Create an instance of the SVG Generator.
      SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
      drawChromatograms(svgGenerator,true,dim,expsToExport,resultsToUse,chromsToUse,analsToExport,analInColumn);
      boolean useCSS = true; // we want to use CSS style attributes
      BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream(fileToStore));
      Writer out = new OutputStreamWriter(stream, "UTF-8");
      svgGenerator.stream(out, useCSS);
    }
  }
  
  private void drawChromatograms(Graphics g, boolean ignorePlatformSpecificSettings, Dimension dim, Vector<String> expsToExport,Vector<String> resultsToUse,Vector<String> chromsToUse, Vector<String> analsToExport, boolean analInColumn){
    // paint the white background
    Graphics2D g2 = (Graphics2D)g;
    g2.setColor(Color.WHITE);
    g2.fillRect(0, 0, dim.width, dim.height);
    int topMargin = 2;
    int marginLegendPicture = 2;
    int bottomMargin = 2;
    int leftMargin = 2;
    int rightMargin = 2;
    double rotationAngle = -Math.PI/2.0;

    // painting the legends
    Vector<String> columnLegend = new Vector<String>(expsToExport);
    Vector<String> rowLegend = new Vector<String>(analsToExport);
    if (analInColumn){
      columnLegend = new Vector<String>(analsToExport);
      rowLegend = new Vector<String>(expsToExport);
    }
    g2.setColor(Color.BLACK);
    Font bigFont = new Font("SansSerif", Font.PLAIN, 200);
    g2.setFont(bigFont);
    FontMetrics bigFontMetrics = g2.getFontMetrics();
    
    Font legendFont = new Font("SansSerif", Font.PLAIN, 12);
    g2.setFont(legendFont);
    FontMetrics legendFontMetrics = g2.getFontMetrics();
    int longestColumnLegend = 0;
    for (String legend : columnLegend){
      int sw = legendFontMetrics.stringWidth(legend);
      if (sw>longestColumnLegend)longestColumnLegend = sw;
    }
    int longestRowLegend = 0;
    for (String legend : rowLegend){
      int sw = legendFontMetrics.stringWidth(legend);
      if (sw>longestRowLegend)longestRowLegend = sw;
    }
    int maxpaintwidth = dim.width-leftMargin-rightMargin-marginLegendPicture-longestRowLegend;
    int maxpaintheight = dim.height-topMargin-bottomMargin-marginLegendPicture-longestColumnLegend;
    int widthOnePicture = maxpaintwidth/columnLegend.size();
    if (widthOnePicture<(legendFontMetrics.getHeight()*8)){
      errorString_ = "The picture width is too small for the selection! Use at least "+(legendFontMetrics.getHeight()*8*columnLegend.size()+leftMargin+rightMargin+marginLegendPicture+longestRowLegend)+" or more!";
      return;
    }
    int heightOnePicture = maxpaintheight/rowLegend.size();
    if (heightOnePicture<(legendFontMetrics.getHeight()*2)){
      errorString_ = "The picture height is too small for the selection! Use at least "+(legendFontMetrics.getHeight()*2*rowLegend.size()+1+topMargin+bottomMargin+marginLegendPicture+longestColumnLegend)+" or more!";
      return;
    }
    int leftPictureMargin = leftMargin+marginLegendPicture+longestRowLegend;
    int topPictureMargin = topMargin+marginLegendPicture+longestColumnLegend;
    for (int i=0; i!=rowLegend.size(); i++){
      String legend = rowLegend.get(i);
      int marginFromTheTop = topMargin+longestColumnLegend+marginLegendPicture;
      int xPos = leftMargin+longestRowLegend-legendFontMetrics.stringWidth(legend);
      int yPos = marginFromTheTop+(i*maxpaintheight)/rowLegend.size()+(heightOnePicture+legendFontMetrics.getHeight())/2;
      g2.drawString(legend, xPos, yPos);
    }
    if (!useMacSpecificSettings(ignorePlatformSpecificSettings))
      g2.rotate(rotationAngle);
    for (int i=0;i!=columnLegend.size();i++){
      String legend = columnLegend.get(i);
      if (useMacSpecificSettings(ignorePlatformSpecificSettings)){       
        BufferedImage image = new BufferedImage(bigFontMetrics.stringWidth(legend), bigFontMetrics.getHeight(),
            BufferedImage.TYPE_INT_ARGB);
        Graphics2D g22 = (Graphics2D) image.getGraphics();
        g22.setColor(Color.BLACK);
        g22.setFont(bigFont);
        g22.drawString(legend,0,bigFontMetrics.getAscent());
        int sw2 = bigFontMetrics.stringWidth(legend);
        int swXExtension2 = (int)(((double)sw2)*Math.cos(-rotationAngle));
        int swYExtension2 = (int)(((double)sw2)*Math.sin(-rotationAngle));
        int shYExtension2 = (int)(((double)bigFontMetrics.getHeight())*Math.cos(-rotationAngle));
        int shXExtension2 = (int)(((double)bigFontMetrics.getHeight())*Math.sin(-rotationAngle));
      
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
        image2 = MacOSResizer.resizeTrick(image2, legendFontMetrics.getHeight(),legendFontMetrics.stringWidth(legend));
        g2.drawImage(image2,(i*maxpaintwidth)/ columnLegend.size()+(widthOnePicture)/2,topMargin+longestColumnLegend-legendFontMetrics.stringWidth(legend),null);

      }else{
        g2.drawString(legend,-topMargin-longestColumnLegend,leftMargin+longestRowLegend+marginLegendPicture + (i*maxpaintwidth)/columnLegend.size()+(widthOnePicture+legendFontMetrics.getHeight())/2);
      }
    }    
    
    if (!useMacSpecificSettings(ignorePlatformSpecificSettings))
      g2.rotate(-rotationAngle);
    
    // now paint the pictures into the boxes
    currentLipidCount_ = 0;
    for (int i=0;i!=chromsToUse.size();i++){
      String chroFile = chromsToUse.get(i);
      String[] filePaths = StringUtils.getChromFilePaths(chroFile);
      try {
        ChromatogramReader reader_ = new ChromatogramReader(filePaths[1], filePaths[2],filePaths[3], filePaths[0],LipidomicsConstants.isSparseData(),LipidomicsConstants.getChromSmoothRange());
        Hashtable<String,Boolean> showMods = new Hashtable<String,Boolean>();
        QuantificationResult result = LDAResultReader.readResultFile(resultsToUse.get(i), showMods);
        Vector<LipidParameterSet> params = result.getIdentifications().get(lipidClass_);
        Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>> resultHash = new Hashtable<String,Hashtable<String,Hashtable<String,LipidParameterSet>>>();
        for (LipidParameterSet set:params){
          String displayString = set.getNameString();
          if (showMods.get(lipidClass_)) displayString+="_"+set.getModificationName();
          String[] molRtAndMod = StaticUtils.extractMoleculeRtAndModFromMoleculeName(displayString);
          if (molRtAndMod[1] == null) molRtAndMod[1] = "";
          if (molRtAndMod[2] == null) molRtAndMod[2] = "";
          Hashtable<String,Hashtable<String,LipidParameterSet>> rtHash = new Hashtable<String,Hashtable<String,LipidParameterSet>>();
          if (resultHash.containsKey(molRtAndMod[0])) rtHash = resultHash.get(molRtAndMod[0]);
          Hashtable<String,LipidParameterSet> modHash = new Hashtable<String,LipidParameterSet>();
          if (rtHash.containsKey(molRtAndMod[1])) modHash = rtHash.get(molRtAndMod[1]);
          modHash.put(molRtAndMod[2], set);
          rtHash.put(molRtAndMod[1], modHash);
          resultHash.put(molRtAndMod[0], rtHash);
        }
        for (int j=0; j!=analsToExport.size(); j++){
          if (interrupt_){
            errorString_ = "The export has been interrupted! The chroms until this point have been exported!";
            this.finished_ = true;
            return;
          }
          String anal = analsToExport.get(j);
          currentLipidCount_++;
          currentLipid_ = anal;
          currentExperiment_ = expsToExport.get(i);
          String[] molRtAndMod = StaticUtils.extractMoleculeRtAndModFromMoleculeName(anal);
          if (molRtAndMod[1] == null) molRtAndMod[1] = "";
          if (molRtAndMod[2] == null) molRtAndMod[2] = "";
          if (resultHash.containsKey(molRtAndMod[0])){
            Hashtable<String,Hashtable<String,LipidParameterSet>> rtHash = resultHash.get(molRtAndMod[0]);
            LipidParameterSet set = null;
            Vector<LipidParameterSet> sets = new Vector<LipidParameterSet>();
            for (String rt : rtHash.keySet()){
              if (rtTolerance_==null || isLookup_.containsKey(anal) || esLookup_.containsKey(anal) || StaticUtils.isWithinTolerance(rtTolerance_,Double.valueOf(rt),Double.valueOf(molRtAndMod[1]))){
                if (rtHash.get(rt).containsKey(molRtAndMod[2])) sets.add(rtHash.get(rt).get(molRtAndMod[2]));
                else if (molRtAndMod[2].equalsIgnoreCase("") && showMods.get(lipidClass_) && rtHash.get(rt).size()==1) sets.add(rtHash.get(rt).values().iterator().next());
              }
            }
            int evidence = StaticUtils.NO_MS2;
            if (sets.size()==1){
              set = sets.get(0);
              evidence = StaticUtils.checkMS2Evidence(set);
            } else if (sets.size()>1) { 
              set = new LipidParameterSet(sets.get(0).Mz[0], sets.get(0).Peptide, sets.get(0).getDoubleBonds(),sets.get(0).getModificationName(),0.0, 
                  sets.get(0).getAnalyteFormula(), sets.get(0).getModificationFormula(), sets.get(0).getCharge(),sets.get(0).getOhNumber());
              set.LowerMzBand = sets.get(0).LowerMzBand;
              set.UpperMzBand = sets.get(0).UpperMzBand;
              Vector<CgProbe> probes = new Vector<CgProbe>();
              for (int l=0; l!=sets.size(); l++){
                LipidParameterSet set0 = sets.get(l);
                int evd0 = StaticUtils.checkMS2Evidence(set0);
                if (evd0>evidence) evidence = evd0; 
                for (int k=0; k!=set0.ProbeCount();k++) probes.add(set0.Probe(k));
              }
              set.setProbes(probes);
            }
            int x0 = leftPictureMargin;
            if (analInColumn)
              x0 += (j*maxpaintwidth)/columnLegend.size();
            else
              x0 += (i*maxpaintwidth)/columnLegend.size();
            int y0 = topPictureMargin;
            if (analInColumn)
              y0 += ((i+1)*maxpaintheight)/rowLegend.size();
            else
              y0 += ((j+1)*maxpaintheight)/rowLegend.size();
            if (set!=null){
              Vector<CgProbe> storedProbes = new Vector<CgProbe>();
              for (int k=0; k!=set.ProbeCount();k++){
                storedProbes.add(set.Probe(k));
              }
              CgChromatogram chrom = reader_.readChromatogram(set.Mz[0]-set.LowerMzBand, set.Mz[0]+set.UpperMzBand,result.getMsLevels().get(lipidClass_));
              chrom.Smooth(LipidomicsConstants.getChromSmoothRange(),
                  LipidomicsConstants.getChromSmoothRepeats());
              chrom.GetMaximumAndAverage();
              float maxValue = chrom.getM_peakValue();
              float maxFound = Analyzer.getHighestIntensityOfProbes(storedProbes,chrom);
              float zoomFactor = 1f;
              if (maxFound>0) zoomFactor = maxValue/maxFound;
              if (evidence>StaticUtils.NO_MS2){
                if (evidence==StaticUtils.PERCENTAL_SPLIT) g.setColor(Color.ORANGE);
                else if (evidence==StaticUtils.PERCENTAL_SPLIT) g.setColor(Color.YELLOW);
                else if (evidence==StaticUtils.MS2_FULL) g.setColor(LipidomicsTableCellRenderer.BRIGHT_GREEN);
                g.fillRect(x0+1, y0-heightOnePicture+1, widthOnePicture-2, heightOnePicture-2);
                g.setColor(Color.BLACK);
              }
              Lipidomics2DPainter.draw2DDiagram(g2, chrom, storedProbes, x0, y0, widthOnePicture, heightOnePicture,
                  MSMapViewer.DISPLAY_TIME_MINUTES,zoomFactor);
            }
            if (LipidomicsConstants.isChromExportShowLegend()){
              String descr = currentExperiment_+" "+currentLipid_;
              int sw = legendFontMetrics.stringWidth(descr);
              g2.drawString(descr, x0+(widthOnePicture-sw)/2,y0-heightOnePicture/2);
            }
          }
        }
      }
      catch (CgException e) {
        File chro = new File(chroFile);
        new WarningMessage(new JFrame(), "Warning", "The file "+chro.getName()+" does not work because of the following reason: "+e.getMessage());
      }catch (ExcelInputFileException e) {
        File chro = new File(chroFile);
        new WarningMessage(new JFrame(), "Warning", "The file "+chro.getName()+" does not work because of the following reason: "+e.getMessage());
      }
    }
  }
  
  private boolean useMacSpecificSettings(boolean ignorePlatformSettings){
    if (!ignorePlatformSettings&&Settings.isOSMacAndJavaLookAndFeel())
      return true;
    else return false;
  }
  
  public void stopThread(){
    this.interrupt_ = true;
  }
  
}
