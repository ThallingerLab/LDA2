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

package at.tugraz.genome.lda;

import java.awt.Button;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.MenuItem;
import java.awt.Panel;
import java.awt.Polygon;
import java.awt.PopupMenu;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.nio.ByteBuffer;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JFrame;

import at.tugraz.genome.dbutilities.Base64;
import at.tugraz.genome.lda.exception.QuantificationException;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.quantification.LipidomicsChromatogram;
import at.tugraz.genome.lda.swing.AreaSettingDialog;
import at.tugraz.genome.lda.swing.RangeColor;
import at.tugraz.genome.lda.vos.AreaSettingVO;
import at.tugraz.genome.maspectras.chromaviewer.MSMapViewer;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.lda.quantification.LipidomicsDefines;
import at.tugraz.genome.maspectras.quantification.CgChromatogram;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.Probe3D;
import at.tugraz.genome.maspectras.utils.Calculator;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class Lipidomics2DPainter extends Panel implements ActionListener,MouseMotionListener,
  MouseListener
{
  private static final long serialVersionUID = 6648556947005652657L;
  private String[] dataStrings_;
  protected Hashtable<Integer,Float> retentionTimes_;
  protected Hashtable<Integer,Float> retentionTimesOriginal_;
  private float mzStart_;
  private float mzStop_;
  private int resolutionFactor_;
  private float stepSize_;
  private ActionListener superiorListener_;
  protected int displayTime_;
  protected boolean raw_;
  private int currentStart_;
  private int currentStop_;
  protected float maxIntensity_;
  protected float timeCorrectionFactor_;
  private float chromMzStart_;
  private float chromMzStop_;
  Vector<CgProbe> storedProbes_;
  Vector<CgProbe> selectedProbes_;
  Vector<CgProbe> removedStoredProbes_;
  
  private Cursor m_crosshair;
  protected TextField m_txtMz;
  protected TextField m_txtTime;
  protected TextField m_txtIntensity;
  protected Button m_plusGain;
  protected Button m_minusGain;
  protected Button m_plusShift;
  protected Button m_minusShift;
  private PopupMenu m_popup;
  private MenuItem m_addUserProbe1;
  private MenuItem m_addUserProbe2;
  private MenuItem m_addUserProbe3;
  private MenuItem m_addUserProbe4;
  private MenuItem m_enterUserProbe;
  private MenuItem m_deleteUserProbe;
  public int lastPopupX = 0;
  public int lastPopupY = 0;
  
  protected float m_2dGain_; // Gain for 2d Display
  protected float m_minDispTime2d_;
  protected float m_maxDispTime2d_;

  protected float m_minTime_, m_maxTime_; // Retention Time boundaries
  
  private int isotopeNumber_;
  
  protected LipidomicsChromatogram cr_;
  private LipidomicsAnalyzer analyzer_;
  private CgProbe probeToAdjust_;
  private boolean isProbeToAdjustSystemProbe_;

  private AreaSettingDialog areaSettings_;
  private int charge_;
  private int msLevel_;
  
  public static final int DISPLAY_TIME_MZ = 100;
  
  protected boolean relativeIntensity_;
  
  private static final float defaultIntensityRangeMultiplicator_ = 1.05f;
  
  private int leftMouseClickLocationX_;
  private int leftMouseClickLocationY_;
  
  private Float ms2Position_;
  
  /** interpolate data points, where MS1 measurements are sparse */
  private boolean doSparseCorrection_ = false;

  public Lipidomics2DPainter(LipidomicsAnalyzer analyzer,String[] dataStrings, Hashtable<Integer,Float> retentionTimes, float mzStart, float mzStop, int resolutionFactor,float stepSize, ActionListener listener,int displayTime,
      float start, float stop, boolean raw, Vector<CgProbe> storedProbes, Vector<CgProbe> selectedProbes,int isotopeNumber, int charge, int msLevel){
    this(analyzer,dataStrings, retentionTimes, null, mzStart, mzStop, resolutionFactor, stepSize, listener, displayTime, start, stop, raw, storedProbes, selectedProbes, isotopeNumber, charge, msLevel);
  }

  
  public Lipidomics2DPainter(LipidomicsAnalyzer analyzer,String[] dataStrings, Hashtable<Integer,Float> retentionTimes,  Hashtable<Integer,Float> rTimesInterpolated, float mzStart, float mzStop, int resolutionFactor,float stepSize, ActionListener listener,int displayTime,
      float start, float stop, boolean raw, Vector<CgProbe> storedProbes, Vector<CgProbe> selectedProbes,int isotopeNumber, int charge, int msLevel){
    super();
    this.analyzer_ = analyzer;
    this.dataStrings_ = dataStrings;
    this.retentionTimes_ = retentionTimes;
    if (LipidomicsConstants.isSparseData() && rTimesInterpolated!=null){
      retentionTimes_ = rTimesInterpolated;
      this.retentionTimesOriginal_ = retentionTimes;
      doSparseCorrection_ = true;
    }
    this.mzStart_ = mzStart;
    this.mzStop_ = mzStop;
    this.resolutionFactor_ = resolutionFactor;
    this.stepSize_ = stepSize;
    this.superiorListener_ = listener;
    if (this.superiorListener_==null){
      this.superiorListener_ = this;
    }
    this.displayTime_ = displayTime;
    timeCorrectionFactor_ = 1f;
    if (this.displayTime_==MSMapViewer.DISPLAY_TIME_MINUTES)
      timeCorrectionFactor_ = 60f;
    this.maxIntensity_ = 0;
    this.m_2dGain_ = 1;
    currentStart_ = Math.round(start*this.resolutionFactor_);
    currentStop_ = Math.round(stop*this.resolutionFactor_);
    m_crosshair = new Cursor(Cursor.CROSSHAIR_CURSOR);
    this.addMouseMotionListener(this);
    this.addMouseListener(this);
    storedProbes_ = storedProbes;
    selectedProbes_ = selectedProbes;
    this.removedStoredProbes_ = new Vector<CgProbe>();
    this.isotopeNumber_ = isotopeNumber;
    charge_ = charge;
    msLevel_ = msLevel;
    relativeIntensity_ = false;
    this.leftMouseClickLocationX_ = -1;
    this.leftMouseClickLocationY_ = -1;
    this.initialize();
  }
  
  public void initialize()
  {
//    this.setLayout(new BorderLayout());
    areaSettings_ = new AreaSettingDialog(new JFrame(), "Probe borders","Enter the borders of the peak",this);
    
    m_txtMz = new TextField();
    m_txtMz.setBackground(Color.WHITE);
    m_txtMz.setEnabled(false);
    m_txtMz.setForeground(Color.BLACK);
    m_txtMz.setVisible(true);
    this.add(m_txtMz);
    
    m_txtTime = new TextField();
    m_txtTime.setBackground(Color.WHITE);
    m_txtTime.setForeground(Color.BLACK);
    m_txtTime.setEnabled(false);
    m_txtTime.setVisible(true);
    this.add(m_txtTime);

    m_txtIntensity = new TextField();
    m_txtIntensity.setBackground(Color.WHITE);
    m_txtIntensity.setForeground(Color.BLACK);
    m_txtIntensity.setEnabled(false);
    m_txtIntensity.setVisible(true);
    this.add(m_txtIntensity);

    m_plusGain = new Button();
    m_plusGain.setLabel("+ Gain");
    m_plusGain.setActionCommand("Plus2dGain");
    m_plusGain.addActionListener(this);
    this.add(m_plusGain);
    m_minusGain = new Button();
    m_minusGain.setLabel("- Gain");
    m_minusGain.setActionCommand("Minus2dGain");
    m_minusGain.addActionListener(this);
    this.add(m_minusGain);

    m_plusShift = new Button();
    m_plusShift.setLabel(">>");
    m_plusShift.setActionCommand("Plus2dShift");
    m_plusShift.addActionListener(this);
    this.add(m_plusShift);
    m_minusShift = new Button();
    m_minusShift.setLabel("<<");
    m_minusShift.setActionCommand("Minus2dShift");
    m_minusShift.addActionListener(this);
    this.add(m_minusShift);
    

    m_popup = new PopupMenu("Peak selection");
    m_addUserProbe4 = new MenuItem("Determine Area (3D)");
    m_addUserProbe4.addActionListener(this.superiorListener_);
    m_popup.add(m_addUserProbe4);
    m_addUserProbe1 = new MenuItem("Determine Area");
    m_addUserProbe1.addActionListener(this.superiorListener_);
    m_popup.add(m_addUserProbe1);
    m_addUserProbe2 = new MenuItem("Determine Area (Col)");
    m_addUserProbe2.addActionListener(this.superiorListener_);
    m_popup.add(m_addUserProbe2);
    m_addUserProbe3 = new MenuItem("Determine Area (Greedy)");
    m_addUserProbe3.addActionListener(this.superiorListener_);
    m_popup.add(m_addUserProbe3);
    m_enterUserProbe = new MenuItem("Enter probe borders");
    m_enterUserProbe.addActionListener(this);
    m_popup.add(m_enterUserProbe);    
    m_popup.addSeparator();
    m_deleteUserProbe = new MenuItem("Delete Area");
    m_deleteUserProbe.addActionListener(this.superiorListener_);
    m_popup.add(m_deleteUserProbe);
    this .add(m_popup);

  ////    hideFields();
    this.setInputLocations();
  }
  
  protected void setInputLocations(){
    m_txtMz.setBounds(this.getWidth() - 100, 5, 95, 20);    
    m_txtTime.setBounds(this.getWidth() - 100, 27, 95, 20);
    m_txtIntensity.setBounds(this.getWidth() - 100, 49, 95, 20);
    m_plusGain.setBounds(this.getWidth() - 65, m_txtIntensity.getY() + 23, 60,
        24);
    m_minusGain.setBounds(this.getWidth() - 65, m_plusGain.getY() + 25, 60, 24);
    m_plusShift.setBounds(this.getWidth() - 35, m_minusGain.getY() + 25, 30, 24);
    m_minusShift.setBounds(this.getWidth() - 65, m_minusGain.getY() + 25, 30,24);

  }

  
  public void preChromatogramExtraxtion(float start, float stop){
    int startAsInt = Math.round(start*this.resolutionFactor_);
    int stopAsInt = Math.round(stop*this.resolutionFactor_);
    boolean setMinMaxTimes = false;
    if (cr_==null)
      setMinMaxTimes = true;
    if (cr_==null || startAsInt != currentStart_||stopAsInt!=currentStop_){
//      System.out.println("Calculate new Chromatogram"+start+";"+stop);
      cr_ = this.extractChromatogram(start,stop);
      m_minTime_ = 0;
      m_maxTime_ = cr_.Value[cr_.Value.length-1][0];
      if (setMinMaxTimes)
        this.m_minDispTime2d_ = m_minTime_;
        this.m_maxDispTime2d_ =  m_maxTime_;
    }else{
//      System.out.println("Use Existing Chromatogram");
    }
  }
  
  public void draw2DDiagram(float start, float stop, boolean raw){
    draw2DDiagram(this.getGraphics(), start, stop, raw);
  }
  
  public void draw2DDiagram(Graphics g){
    draw2DDiagram(g,Float.valueOf(String.valueOf(currentStart_))/this.resolutionFactor_, Float.valueOf(String.valueOf(currentStop_))/this.resolutionFactor_, this.raw_);
  }
  
  protected void draw2DDiagram(Graphics g, float start, float stop, boolean raw){
    chromMzStart_ = start;
    chromMzStop_ = stop;
    raw_ = raw;
    int x0 = leftMargin_2d();
    int y0 = this.getHeight() - bottomMargin_2d();
    int w0 = this.diagramWidth_2d();
    int h0 = this.diagramHeight_2d();
    int startAsInt = Math.round(start*this.resolutionFactor_);
    int stopAsInt = Math.round(stop*this.resolutionFactor_);
    boolean setMinMaxTimes = false;
    if (cr_==null)
      setMinMaxTimes = true;
    if (cr_==null || startAsInt != currentStart_||stopAsInt!=currentStop_){
//      System.out.println("Calculate new Chromatogram");
      cr_ = this.extractChromatogram(start,stop);
      m_minTime_ = 0;
      m_maxTime_ = cr_.Value[cr_.Value.length-1][0];
      if (setMinMaxTimes)
        this.m_minDispTime2d_ = m_minTime_;
        this.m_maxDispTime2d_ =  m_maxTime_;
    }else{
//      System.out.println("Use Existing Chromatogram");
    }
    Vector<CgProbe> paintableProbes = new Vector<CgProbe>();
    for (int i = 0; i != this.storedProbes_.size(); i++) {
      CgProbe aProbe = (CgProbe) this.storedProbes_.get(i);
      boolean removed = false;
      for (int j = 0; j != this.removedStoredProbes_.size(); j++) {
        if (this.sameCgProbe((CgProbe) this.removedStoredProbes_.get(j), aProbe)) {
          removed = true;
        }
      }
      if (!removed)
        paintableProbes.add(aProbe);
    }


    paint2dAreas(paintableProbes, Color.RED);
    paint2dAreas(this.selectedProbes_, Color.GREEN);

    draw2DDiagram(g, cr_, x0,y0,w0,h0, m_minDispTime2d_, m_maxDispTime2d_,maxIntensity_,m_2dGain_,
        displayTime_, timeCorrectionFactor_, raw_,true,relativeIntensity_,false,this.ms2Position_,false);
    currentStart_ = startAsInt;
    currentStop_ = stopAsInt;
    this.m_txtMz.setText("mz="+Calculator.roundFloat((stop+start)/2f, 3)+"D");

  }
  
  public static void draw2DDiagram(Graphics gx, CgChromatogram cr, Vector<CgProbe> probes, int xIn, int yIn, int wIn, int hIn,
      int displayTime, float zoomFactor){
    float maxIntensity = 0;
    for (int scanNumber=0; scanNumber!=cr.Value.length; scanNumber++){
      if (cr.Value[scanNumber][1]> maxIntensity){
        maxIntensity = cr.Value[scanNumber][1];
      }
    }
    float timeCorrectionFactor = 1f;
    if (displayTime==MSMapViewer.DISPLAY_TIME_MINUTES)
      timeCorrectionFactor = 60f;
    int x0 = xIn+4;
    int y0 = yIn-7;
    int w0 = wIn-15;
    int h0 = hIn-8;
    
    paint2dAreas(gx,cr, probes, Color.RED, 0, x0, y0, w0, h0, 0,cr.Value[cr.Value.length-1][0],maxIntensity,zoomFactor,false,false);
    draw2DDiagram(gx, cr, x0, y0, w0, h0, 0, cr.Value[cr.Value.length-1][0], maxIntensity, zoomFactor, displayTime, timeCorrectionFactor,
        false,false,false,false,null,true);
  }

  /**
   * draws a 2D diagram with the provided Graphics component gx
   * @param gx Graphics component for drawing
   * @param cr chromatogram to be drawn
   * @param x0 diagram start position (x)
   * @param y0 diagram start position (y)
   * @param w0 maximum width of diagram in pixels
   * @param h0 maximum height of diagram in pixels
   * @param m_minDispTime2d start of range of diagram (in its original units)
   * @param m_maxDispTime2d stop of range of diagram (in its original units)
   * @param maxIntensity highest intensity available (for drawing of coordinates)
   * @param m_2dGain zoom factor
   * @param displayTime x axis type (MSMapViewer.DISPLAY_TIME_SECONDS, MSMapViewer.DISPLAY_TIME_MINUTES, DISPLAY_TIME_MZ)
   * @param timeCorrectionFactor for correcting the time value
   * @param raw shall the raw intensity values be taken for display or the smoothed ones
   * @param drawLegend draw a legend in the diagram
   * @param relativeValue should relative or absolute values be used
   * @param barChart draw a bar chart or a xy line chart
   * @param ms2Position paint a line where the currently displayed MS2 spectrum is in the chromatogram
   * @param stopAtH0 for chromatogram exporting = true - if false the values may be painted into another chromatogram 
   */
  protected static void draw2DDiagram(Graphics gx, CgChromatogram cr, int x0, int y0, int w0, int h0,
      float m_minDispTime2d, float m_maxDispTime2d, float maxIntensity, float m_2dGain, int displayTime, float timeCorrectionFactor, boolean raw,
      boolean drawLegend, boolean relativeValue, boolean barChart, Float ms2Position, boolean stopAtH0){
    draw2DDiagram(gx, cr, x0, y0, w0, h0, m_minDispTime2d, m_maxDispTime2d, maxIntensity, m_2dGain, displayTime, timeCorrectionFactor, raw,
        drawLegend, relativeValue, barChart, ms2Position, null, stopAtH0);
  }

  
  /**
   * draws a 2D diagram with the provided Graphics component gx
   * @param gx gx Graphics component for drawing
   * @param cr chromatogram to be drawn
   * @param x0 diagram start position (x)
   * @param y0 diagram start position (y)
   * @param w0 maximum width of diagram in pixels
   * @param h0 maximum height of diagram in pixels
   * @param m_minDispTime2d start of range of diagram (in its original units)
   * @param m_maxDispTime2d stop of range of diagram (in its original units)
   * @param maxIntensity highest intensity available (for drawing of coordinates)
   * @param m_2dGain zoom factor
   * @param displayTime x axis type (MSMapViewer.DISPLAY_TIME_SECONDS, MSMapViewer.DISPLAY_TIME_MINUTES, DISPLAY_TIME_MZ)
   * @param timeCorrectionFactor for correcting the time value
   * @param raw shall the raw intensity values be taken for display or the smoothed ones
   * @param drawLegend draw a legend in the diagram
   * @param relativeValue should relative or absolute values be used
   * @param barChart draw a bar chart or a xy line chart
   * @param ms2Position paint a line where the currently displayed MS2 spectrum is in the chromatogram
   * @param rangeColors shall certain ranges of the chromatogram be painted in a different color
    * @param stopAtH0 for chromatogram exporting = true - if false the values may be painted into another chromatogram 
  */
  protected static void draw2DDiagram(Graphics gx, CgChromatogram cr, int x0, int y0, int w0, int h0,
      float m_minDispTime2d, float m_maxDispTime2d, float maxIntensity, float m_2dGain, int displayTime, float timeCorrectionFactor, boolean raw,
      boolean drawLegend, boolean relativeValue, boolean barChart, Float ms2Position, Vector<RangeColor> rangeColors, boolean stopAtH0)
  {

    int x = 0, y = 0;
    int i, j;
//    boolean firstPoint;
    String s;
    float maxInt = getMaxInt(maxIntensity,relativeValue);
    // **** Draw the Areas ****
    
//    len = cr.Value.length;

    // **** Draw the diagram ****
    drawDiagramItself(gx, cr, x0, y0, w0, h0, m_minDispTime2d,m_maxDispTime2d, maxIntensity, m_2dGain, raw, relativeValue, barChart, rangeColors, stopAtH0);
    gx.setColor(Color.BLACK);

    // **** Draw the x coordinates ****

    float t = 0;
    float dt;

    gx.setFont(new Font("SansSerif", Font.PLAIN, 10));
    FontMetrics fm = gx.getFontMetrics();
    int sw;

    gx.drawLine(x0, y0, x0 + w0 + 3, y0);
    gx.drawLine(x0 + w0 + 3, y0 - 2, x0 + w0 + 3, y0 + 2);
    gx.drawLine(x0 + w0 + 3, y0 - 2, x0 + w0 + 9, y0);
    gx.drawLine(x0 + w0 + 3, y0 + 2, x0 + w0 + 9, y0);

    gx.drawLine(x0, y0, x0, y0 - h0 + 5);
    gx.drawLine(x0 - 2, y0 - h0 + 5, x0 + 2, y0 - h0 + 5);
    gx.drawLine(x0 - 2, y0 - h0 + 5, x0, y0 - h0 - 1);
    gx.drawLine(x0 + 2, y0 - h0 + 5, x0, y0 - h0 - 1);

    dt = 100;
    t = (m_maxDispTime2d - m_minDispTime2d) / 10;
    t = t/timeCorrectionFactor;
    if (t < 1)
      dt = 1;
    else if (t < 10)
      dt = 10;
    else if (t < 20)
      dt = 20;
    else if (t < 50)
      dt = 50;
    else if (t < 100)
      dt = 100;
    else if (t < 200)
      dt = 200;
    else if (t < 500)
      dt = 500;
    else if (t < 1000)
      dt = 1000;
    else if (t < 2000)
      dt = 2000;
    else if (t < 5000)
      dt = 5000;
    else if (t < 10000)
      dt = 10000;
    else if (t < 20000)
      dt = 20000;
    else if (t < 50000)
      dt = 50000;
    else
      dt = 100000;
    
    String timeUnit = "Time / s";
    if (displayTime==MSMapViewer.DISPLAY_TIME_MINUTES)
      timeUnit = "Time / min";
    else if (displayTime==DISPLAY_TIME_MZ)
      timeUnit = "m/z";
    j = fm.stringWidth(timeUnit);
    
    i = (int) ((m_minDispTime2d/timeCorrectionFactor) / dt);
    if (i > 0){
      t = (i - 1) * dt;
    }else{
      t = 0;
    }
    while ((t*timeCorrectionFactor) < m_maxDispTime2d) {
      if (t*timeCorrectionFactor >= m_minDispTime2d) {
        x = x0
            + (int) ((w0 * (t* timeCorrectionFactor - m_minDispTime2d))  / (m_maxDispTime2d - m_minDispTime2d));
        y = y0;
        gx.drawLine(x, y, x, y + 3);
        s = Integer.toString((int) t);
        sw = fm.stringWidth(s);

        if ((x + sw / 2 + 10) < (x0 + w0 - j)){
          if (drawLegend) gx.drawString(s, x - sw / 2, y + 12);
        }  
      }
      t += dt;
    }
    if (drawLegend) gx.drawString(timeUnit, x0 + w0 - j + 9, y + 12);

    // **** Draw the y coordinates ****

    float in = 0;
    float din = 0;

    in = maxInt  / m_2dGain / 5;
    if (in < 0.1f)
      din = 0.1f;
    else if (in < .25f)
      din = 0.25f;
    else if (in < 0.5f)
      din = 0.5f;
    else if (in < 1)
      din = 1;
    else if (in < 2.5f)
      din = 2.5f;
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
    else
      din = 1000000000;

    in = 0;
    y = y0 - (int) (h0 * in / maxInt  * m_2dGain);
    // while(in<m_maxIntensity)
    while (y > (y0 - h0)) {
      x = x0;
      y = y0 - (int) (h0 * in / maxInt   * m_2dGain);
      if (y < (y0 - h0 + 10))
        break;
      if (in > 2147483000)
        break;

      gx.drawLine(x, y, x - 3, y);
      s = Integer.toString((int) in);
      if (in<1&&in>0||(in>1&&in<2)||(in>2&&in<3)||(in>7&&in<8)) s = Float.toString(in);
      sw = fm.stringWidth(s);
      if (drawLegend){
        if (in == 0)
          gx.drawString(s, x - sw - 4, y);
        else
          gx.drawString(s, x - sw - 4, y + 4);
      }
      in += din;
    }
    if (relativeValue){
      j = fm.stringWidth("Relative abundance [%]");
      double standRotAngle = -Math.PI/2;
      ((Graphics2D)gx).rotate(standRotAngle);
      gx.drawString("Relative abundance [%]",(h0-j)/2-y0, x/2+fm.getHeight()/2);
      ((Graphics2D)gx).rotate(-standRotAngle);
    }else{
      j = fm.stringWidth("Arb.Units / AU");
      if (drawLegend) gx.drawString("Arb.Units / AU", x - j - 4, y0 - h0 + 6);      
    }
    
    if (ms2Position !=null){
      float values[][] = new float[1][4];
      values[0][0] = ms2Position;
      if (displayTime==MSMapViewer.DISPLAY_TIME_MINUTES) values[0][0] = ms2Position*60f;
      values[0][1] = maxIntensity;
      values[0][2] = maxIntensity;
      int [] coords = getCoordinatesInDiagram(values, 0, x0, y0, w0, h0, m_minDispTime2d, m_maxDispTime2d, maxIntensity, m_2dGain, raw, relativeValue,stopAtH0);
      gx.setColor(Color.BLUE);
      gx.drawLine(coords[0], y0, coords[0], coords[1]);
      gx.setColor(Color.BLACK);
    }
//    m.setMz(cr.getMz());
  }
  
  private LipidomicsChromatogram extractChromatogram(float start, float stop){
    int startIndex = this.getDataIndex(start);
    if (startIndex<1)
      startIndex = 0;
    int stopIndex = this.getDataIndex(stop);
    if (stopIndex>=this.dataStrings_.length)
      stopIndex = this.dataStrings_.length;
    Hashtable<Integer,Float> rts = this.retentionTimes_;
    if (doSparseCorrection_) rts = this.retentionTimesOriginal_;
    LipidomicsChromatogram chrom = new LipidomicsChromatogram(rts.size());
    chrom.LowerMzBand = this.stepSize_/2;
    chrom.UpperMzBand = this.stepSize_/2;
    chrom.Mz = (start+stop)/2;
    for (int i=0;i!=rts.size();i++){
      chrom.Value[i][0] = (rts.get(new Integer(i))).floatValue();
      chrom.Value[i][1] = 0;
    }
    if (dataStrings_!=null){
      for (int i=startIndex; i!=stopIndex;i++){
        if (dataStrings_[i]!=null&&dataStrings_[i].length()>0){
          ByteBuffer buffer = ByteBuffer.wrap(Base64.decode(dataStrings_[i]));
          while (buffer.hasRemaining()){
            int scanNumber = buffer.getInt();
            float intensity = buffer.getFloat();
            chrom.Value[scanNumber][1] += intensity;
            if (chrom.Value[scanNumber][1]> this.maxIntensity_){
              maxIntensity_ = chrom.Value[scanNumber][1];
            }
          }
        }
      }
    }
//    if (sparseData) chrom = new LipidomicsChromatogram(ChromatogramReader.interpolateChromValues(chrom, this.retentionTimes_));
    if (doSparseCorrection_) chrom.doChromValueInterpolation(retentionTimes_);
    chrom.Smooth(LipidomicsConstants.getChromSmoothRange(),
        LipidomicsConstants.getChromSmoothRepeats());
    chrom.GetMaximumAndAverage();
    return chrom;

  }
  
  protected static void drawDiagramItself(Graphics gx, CgChromatogram cr, int x0, int y0, int w0, int h0, float m_minDispTime2d, 
      float m_maxDispTime2d, float maxIntensity, float m_2dGain, boolean raw, boolean relativeValue, boolean barChart, Vector<RangeColor> rangeColors,
      boolean stopAtH0){
    int x = 0, y = 0;
    int xx = 0, yy = 0;
    gx.setColor(Color.BLACK);
    boolean firstPoint = true;
    int len = cr.Value.length;
    for (int i = 0; i < len; i++) {

      if (cr.Value[i][0] < m_minDispTime2d)
        continue;
      if (cr.Value[i][0] > m_maxDispTime2d)
        continue;
      int[] coords = getCoordinatesInDiagram(cr.Value, i, x0, y0, w0, h0, m_minDispTime2d, 
          m_maxDispTime2d, maxIntensity, m_2dGain, raw, relativeValue,stopAtH0);
      
      x = coords[0];
      y = coords[1];
      if (firstPoint == true) {
        xx = x;
        yy = y;
        firstPoint = false;
      }
      if (rangeColors!=null && rangeColors.size()>0){
        RangeColor rangeVO = getColorAccordingToXPos(cr.Value[i][0],rangeColors);
        if (rangeVO!=null){
          gx.setColor(rangeVO.getColor());
        }else gx.setColor(Color.BLACK);
      }
      if (barChart){
        gx.drawLine(x, y0, x, y<0 ? 0 : y);
      }else
        gx.drawLine(xx, yy, x, y);
      xx = x;
      yy = y;
    }
    gx.setColor(Color.BLACK);
  }
  
  /**
   * checks if the position (time) is inside a found hit and returns the corresponding color
   * @param x the retention time position
   * @param rangeColors the ranges of found hits and their colors
   * @return the found RangeColorVO or null if nothing is found
   */
  protected static RangeColor getColorAccordingToXPos(float x,Vector<RangeColor> rangeColors){
    for (RangeColor vo : rangeColors){
      if (vo.insideRange(x)) return vo;
    }
    return null;
  }
    
  public void paint(Graphics g)
  {
    draw2DDiagram(g);
    setInputLocations();
    showFields();
  }

  @SuppressWarnings("deprecation")
  public void reshape(int x, int y, int width, int height){
    this.hideFields();
    super.reshape(x, y, width, height);
  }  
  public void hideFields()
  {
    if (m_txtTime.isVisible()){
      if (m_txtMz!=null) m_txtMz.setVisible(false);
      m_txtTime.setVisible(false);
      m_txtIntensity.setVisible(false);
      m_plusGain.setVisible(false);
      m_minusGain.setVisible(false);
      m_plusShift.setVisible(false);
      m_minusShift.setVisible(false);
    }
  }
  
  public void showFields()
  {
    if (!m_txtTime.isVisible()){
      if (m_txtMz!=null) m_txtMz.setVisible(true);
      m_txtTime.setVisible(true);
      m_txtIntensity.setVisible(true);
      m_plusGain.setVisible(true);
      m_minusGain.setVisible(true);
      m_plusShift.setVisible(true);
      m_minusShift.setVisible(true);
    }
  }
  
  private int getDataIndex(float mzValue){
    String toCut = String.valueOf(Calculator.roundFloat((mzValue - this.mzStart_)*this.resolutionFactor_,0));
    return Integer.parseInt(toCut.substring(0,toCut.indexOf(".")));
  }
  
  // =============================================================================================
  // Helpers for the 2d Panel:
  // =============================================================================================

  protected int leftMargin_2d()
  {
    return 75;
  }

  private int rightMargin_2d()
  {
    return 15;
  }

  private static int topMargin_2d()
  {
    return 5;
  }

  protected int bottomMargin_2d()
  {
    return 20;
  }

  public void actionPerformed(ActionEvent e)
  {
    String command = e.getActionCommand();
    if (command.equalsIgnoreCase("Plus2dGain")){
      this.m_2dGain_ = this.m_2dGain_ * 1.1f;
      this.update(this.getGraphics());
//      this.paint(this.getGraphics());
    }
    if (command.equalsIgnoreCase("Minus2dGain")){
      this.m_2dGain_ = this.m_2dGain_/1.1f;
      this.update(this.getGraphics());      
    }
    if (command.equalsIgnoreCase("Plus2dShift")){
      this.set2dShift(+1);
      this.update(this.getGraphics());      
    }
    if (command.equalsIgnoreCase("Minus2dShift")){
      this.set2dShift(-1);
      this.update(this.getGraphics());      
    }
    if (e.getActionCommand() == "Determine Area") {
//    removeOverlappingProbes(cx);
//    if (restoreSystemArea(cx) == false)
//      determineArea(m_2dPnl.lastPopupX, 0);
//      determineArea(this.lastPopupX, 0);
      this.determineArea(LipidomicsDefines.StandardValleyMethod);
    return;
  }else if (e.getActionCommand() == "Determine Area (Col)") {
    this.determineArea(LipidomicsDefines.EnhancedValleyMethod);
    return;
  }else if (e.getActionCommand() == "Determine Area (Greedy)") {
    this.determineArea(LipidomicsDefines.GreedySteepnessReductionMethod);
    return;
  } else if (e.getActionCommand() == "Delete Area") {
    this.deleteArea();
  } else if (e.getActionCommand().equalsIgnoreCase("Enter probe borders")) {
    CgProbe usrPr = getUserProbe(lastPopupX);
    CgProbe sysPr = getSystemProbe(lastPopupX);
    boolean threeD = false;
    float mzStart = chromMzStart_;
    float mzStop = chromMzStop_;
    float time = getRetentionTime(lastPopupX);
    float timeStart = time-60;
    float timeStop = time+60;
    if (usrPr!=null || sysPr!=null){
      if (usrPr!=null){
        probeToAdjust_ = usrPr;
        isProbeToAdjustSystemProbe_ = false;
      }else{
        probeToAdjust_ = sysPr;
        isProbeToAdjustSystemProbe_ = true;
      }
      timeStart = probeToAdjust_.LowerValley;
      timeStop = probeToAdjust_.UpperValley;
      if (probeToAdjust_ instanceof Probe3D){
        threeD = true;
        Probe3D probe3D = (Probe3D)probeToAdjust_;
        mzStart = probe3D.getEllipseMzPosition()-probe3D.getEllipseMzStretch();
        mzStop = probe3D.getEllipseMzPosition()+probe3D.getEllipseMzStretch();
        timeStart = probe3D.getEllipseTimePosition()-probe3D.getEllipseTimeStretch();
        timeStop = probe3D.getEllipseTimePosition()+probe3D.getEllipseTimeStretch();
      }
    } else
      probeToAdjust_ = null;
    if (timeStart<m_minDispTime2d_) timeStart = m_minDispTime2d_;
    if (timeStop>m_maxDispTime2d_) timeStop = m_maxDispTime2d_;
    timeStart /= timeCorrectionFactor_;
    timeStop /= timeCorrectionFactor_;
    areaSettings_.setInputFields(threeD, timeStart, timeStop, mzStart, mzStop);
    areaSettings_.setVisible(true);
  }else if (e.getActionCommand().equalsIgnoreCase("acceptAreaSettings")){
    this.acceptAreaSettings();
    if (this.superiorListener_!=null) this.superiorListener_.actionPerformed(e);
  }else if (e.getActionCommand().equalsIgnoreCase("discardAreaSettings")){
    areaSettings_.setVisible(false);
  }  

    
  }
  
  private void acceptAreaSettings(){
    AreaSettingVO settings = areaSettings_.getInputValues();
    if ((settings.getStartTime()*timeCorrectionFactor_)<0)
      new WarningMessage(new JFrame(), "Error", "Your start-time must be greater than 0!");
    else if ((settings.getStopTime()*timeCorrectionFactor_)>(m_maxDispTime2d_+30))
      new WarningMessage(new JFrame(), "Error", "Your stop-time must not be greater than "+m_maxDispTime2d_/timeCorrectionFactor_+"! This is outside the range");
    else if (settings.getStartTime()>=settings.getStopTime())
      new WarningMessage(new JFrame(), "Error", "Your start-time must not be greater or equal your stop-time!");
    else if (settings.getStartMz()<mzStart_)
      new WarningMessage(new JFrame(), "Error", "Your start-m/z must not be lower than "+mzStart_+"! This is outside the range");
    else if (settings.getStopMz()>mzStop_)
      new WarningMessage(new JFrame(), "Error", "Your stop-m/z must not be greater than "+mzStop_+"! This is outside the range");
    else if (settings.getStartMz()>=settings.getStopMz())
      new WarningMessage(new JFrame(), "Error", "Your start-m/z must not be greater or equal your stop-m/z!");
    else {
      CgProbe probeToAdd = this.createManuallyAddedProbe(settings);
      areaSettings_.setVisible(false);
      probeToAdd.isotopeNumber = isotopeNumber_;
      addProbeToTheCache(probeToAdd);
    }
    
  }

  public void setRaw(boolean raw)
  {
    this.raw_ = raw;
  }
  
  public boolean setMinDispTime2d(float newVal)
  {
    float newVal2 = newVal*this.timeCorrectionFactor_;
    if (newVal2 < 0)
      return false;
    if (newVal2 < m_minTime_)
      newVal2 = m_minTime_;
    m_minDispTime2d_ = newVal2;
    if (m_minDispTime2d_ < m_maxTime_ / 100)
      m_minDispTime2d_ = 0;
    return true;
  }
  
  public boolean setMaxDispTime2d(float newVal)
  {
    float newVal2 = newVal*this.timeCorrectionFactor_;
    if (newVal2 < m_minDispTime2d_)
      return false;
    if (newVal2 > m_maxTime_)
      newVal2 = m_maxTime_;
    m_maxDispTime2d_ = newVal2;
    return true;
  }

  public void zoomAll()
  {
    m_minDispTime2d_ = m_minTime_;
    m_maxDispTime2d_ = m_maxTime_;
    if (m_minDispTime2d_ < m_maxTime_ / 100)
      m_minDispTime2d_ = 0;
  }
  
  public void mouseDragged(MouseEvent e){}

  public void mouseMoved(MouseEvent e){
    this.setShownIntensity(getIntensity(e.getY()));
    this.setShownRetentionTime(getRetentionTime(e.getX()));
    if (leftMouseClickLocationX_>-1){
      
    }
  }

  public void mouseClicked(MouseEvent e){}

  public void mouseEntered(MouseEvent e)
  {
    setCursor(m_crosshair);
  }

  public void mouseExited(MouseEvent e)
  {
    setCursor(Cursor.getDefaultCursor());
    m_txtIntensity.setText("");
    m_txtTime.setText("");
  }
  
  public void mousePressed(MouseEvent e){
    if (e.getButton()==MouseEvent.BUTTON1){
      leftMouseClickLocationX_ = e.getX();
      leftMouseClickLocationY_ = e.getY();
      paintTriangle(leftMouseClickLocationX_,leftMouseClickLocationY_);
    }
  }
  
  private void paintTriangle(int x, int y){
    Graphics graphics = this.getGraphics();
    triangleWrite(graphics, x, y);  }
  
  private void eraseTriangle(int x, int y){
    Graphics graphics = this.getGraphics();
    graphics.setColor(Color.white);
    triangleWrite(graphics, x, y);
    graphics.setColor(Color.black);
  }
  
  private void triangleWrite(Graphics graphics, int x, int y){
    graphics.drawLine(leftMouseClickLocationX_-8, leftMouseClickLocationY_-5, leftMouseClickLocationX_, leftMouseClickLocationY_);
    graphics.drawLine(leftMouseClickLocationX_+8, leftMouseClickLocationY_-5, leftMouseClickLocationX_, leftMouseClickLocationY_);
    graphics.drawLine(leftMouseClickLocationX_-8, leftMouseClickLocationY_-5, leftMouseClickLocationX_+8, leftMouseClickLocationY_-5);    
  }
  

  public void mouseReleased(MouseEvent e)
  {
    if (leftMouseClickLocationX_>-1){
      eraseTriangle(leftMouseClickLocationX_,leftMouseClickLocationY_);
      leftMouseClickLocationX_ = -1;
      leftMouseClickLocationY_ = -1;
    }
    if (e.isPopupTrigger()||e.getButton()==MouseEvent.BUTTON3) {
      if (this.checkPopupClick(e.getX(), e.getY()) == 0)
        return;
      boolean usrPr = this.anyUserProbeThere(e.getX());
      boolean sysPr = this.anySystemProbeThere(e.getX());
      m_deleteUserProbe.setEnabled(usrPr || sysPr);
      m_addUserProbe1.setEnabled(sysPr == false && usrPr == false);
      m_addUserProbe2.setEnabled(sysPr == false && usrPr == false);
      m_addUserProbe3.setEnabled(sysPr == false && usrPr == false);
      m_addUserProbe4.setEnabled(sysPr == false && usrPr == false);

      this.setCursor(Cursor.getDefaultCursor());
      lastPopupX = e.getX();
      lastPopupY = e.getY();
      m_popup.show(e.getComponent(), e.getX(), e.getY());
    }
  }
  
  protected void setShownIntensity(float newVal)
  {
    if (newVal >= 0) {
      float value = newVal;
      value = getCorrespondingValue(value,maxIntensity_,raw_,relativeIntensity_);
      String s = String.format("Int = %4.3e", value);
      s = s.replace("+", "").replace("e0", "e");
      m_txtIntensity.setText(s);
    } else
      m_txtIntensity.setText("");
  }

  protected float getIntensity(int y)
  {
    if (y < topMargin_2d())
      return -1;
    int h0 = this.diagramHeight_2d();
    if (y > (topMargin_2d() + h0))
      return -1;
    return (h0 + topMargin_2d() - y) * maxIntensity_ * defaultIntensityRangeMultiplicator_ / m_2dGain_ / h0;
  }
  
  protected float getRetentionTime(int x)
  {
    if (x < leftMargin_2d())
      return -1;
    int w0 = this.diagramWidth_2d();
    if (x > leftMargin_2d() + w0)
      return -1;
    return m_minDispTime2d_+ ((x - leftMargin_2d()) * (m_maxDispTime2d_ - m_minDispTime2d_))
        / (float) w0;
  }
  
  public void setShownRetentionTime(float newVal)
  {
    if (newVal >= 0) {
      float time = newVal/this.timeCorrectionFactor_;
      String value = "";
      if (this.displayTime_==MSMapViewer.DISPLAY_TIME_SECONDS)
        value = String.valueOf(Calculator.roundFloat(time, 0));
      if (this.displayTime_==MSMapViewer.DISPLAY_TIME_MINUTES)
        value = String.valueOf(Calculator.roundFloat(time, 2));
      
      String s = value;
      m_txtTime.setText("t = " + s);
    } else
      m_txtTime.setText("");
  }
  
  private void set2dShift(int direction)
  {
    float x, diff, delta;
    diff = m_maxDispTime2d_ - m_minDispTime2d_;
    delta = diff / 40;
    if (direction < 0) {
      x = m_minDispTime2d_ - delta;
      if (x < 0)
        x = 0;
      m_minDispTime2d_ = x;
      m_maxDispTime2d_ = x + diff;
    } else {
      x = m_maxDispTime2d_ + delta;
      if (x > m_maxTime_)
        x = m_maxTime_;
      m_maxDispTime2d_ = x;
      m_minDispTime2d_ = x - diff;
    }
  }
  
  public void nextChromatogram(){
    int diff = this.currentStop_- this.currentStart_;
    if (Float.valueOf((this.currentStop_+diff))/this.resolutionFactor_<=this.mzStop_){
      this.currentStart_ = this.currentStart_+diff;
      this.currentStop_ = this.currentStop_+diff;
      this.cr_ = null;
    }  
  }
  
  public void previousChromatogram(){
    int diff = this.currentStop_- this.currentStart_;
    if (Float.valueOf((this.currentStart_-diff))/this.resolutionFactor_>=this.mzStart_){
        this.currentStart_ = this.currentStart_-diff;
        this.currentStop_ = this.currentStop_-diff;
        this.cr_ = null;
    }
  }
  
  public float getLowerMz(){
    return (Float.valueOf(this.currentStart_)/this.resolutionFactor_);
  }
  
  public float getUpperMz(){
    return (Float.valueOf(this.currentStop_)/this.resolutionFactor_);
  }
  
  private void paint2dAreas(Vector<CgProbe> areas, Color col){
    int x0 = leftMargin_2d();
    int w0 = diagramWidth_2d();
    int y0 = topMargin_2d() + diagramHeight_2d();
    int h0 = diagramHeight_2d();
    paint2dAreas(this.getGraphics(),cr_, areas, col, isotopeNumber_, x0, y0, w0, h0, m_minDispTime2d_,m_maxDispTime2d_,maxIntensity_,m_2dGain_,raw_,true);
  }

  private static void paint2dAreas(Graphics gx, CgChromatogram cr, Vector<CgProbe> areas, Color col, int isotopeNumber,int x0, int y0, int w0, int h0,
      float m_minDispTime2d, float m_maxDispTime2d, float maxIntensity, float m_2dGain, boolean raw, boolean displayLegend)
  {
    int j, len;
    int x = 0, y = 0;
    boolean firstPoint;
    Polygon p;
    String s;
    
    float maxInt = maxIntensity;
    maxInt = maxInt*defaultIntensityRangeMultiplicator_;

    len = cr.Value.length;
    for (CgProbe cp: areas) {
      if (cp == null)
        continue;
      if (cp.AreaStatus != CgAreaStatus.OK)
        continue;
//      System.out.println("IsotopeNumber: "+cp.isotopeNumber);
////      if (cp.Mz < (cr_.Mz + this.getThreshold())
////          && cp.Mz > (cr_.Mz - this.getThreshold())) {
      if (cp.isotopeNumber==isotopeNumber){
        // Paint the Area of the Peak accordingly:
        p = new Polygon();
        firstPoint = true;
        int nearestScanNumber = 0;
        for (j = 0; j < len; j++) {
          if (cr.Value[j][0] < m_minDispTime2d)
            continue;
          if (cr.Value[j][0] < cp.LowerValley)
            continue;
          if (cr.Value[j][0] > m_maxDispTime2d)
            continue;
          if (cr.Value[j][0] > cp.UpperValley)
            continue;
          if (cr.Value[j][0]<cp.Peak)
            nearestScanNumber = j;
          x = x0
              + (int) (w0 * (cr.Value[j][0] - m_minDispTime2d) / (m_maxDispTime2d - m_minDispTime2d));

          float value = getCorrespondingValue(cr.Value[j],maxIntensity,raw,false);
          y = y0 - (int) (h0 * value * m_2dGain / maxInt);

          if (firstPoint == true) {
            p.addPoint(x, y0);
            firstPoint = false;
          }
          p.addPoint(x, y);
        }
        if (p.npoints > 0) {
          p.addPoint(x, y0);
          gx.setColor(col);
          gx.fillPolygon(p);

          // draw the Area:

          x = x0
              + (int) (w0 * (cp.Peak - m_minDispTime2d) / (m_maxDispTime2d - m_minDispTime2d));
          float value = getCorrespondingValue(cr.Value[nearestScanNumber],maxIntensity,raw,false);
          y = y0 - (int) (h0 * value * m_2dGain / maxInt) - 5;

          if (y < 30) {
            y = 30;
            x += 10;
          }

          s = String.format("A=%4.3e", cp.Area);
          s = s.replace("+", "")+" @ "+Calculator.roundFloat(cp.LowerValley/60f,2)+"<t<"+Calculator.roundFloat(cp.UpperValley/60f,2);
              //+ String.format(" @ %.0f<t<%.0f", Calculator.roundFloat(cp.LowerValley/60f,2), Calculator.roundFloat(cp.UpperValley/60f,2));
//          s += String.format(" / -%.2fd+%.2fd", cp.LowerMzBand, cp.UpperMzBand);
//          s += String.format(" / Bckgnd/A=%.1f%%", 100 * cp.Background
//              / cp.Area);
          s = s.replace("e0", "e").replace(",", ".");
          gx.setColor(Color.BLUE);
          if (displayLegend) gx.drawString(s, x, y);
        }
      }
    }
  }

  protected int diagramWidth_2d()
  {
    return this.getWidth() - leftMargin_2d() - rightMargin_2d();
  }
  
  protected int diagramHeight_2d()
  {
    return this.getHeight() - topMargin_2d() - bottomMargin_2d();
  }
  
  private double getThreshold()
  {
    return (this.stepSize_/50);
  }
  
  private boolean in2dDiagram(int x, int y)
  {
    if (x < leftMargin_2d())
      return false;
    if (x > (leftMargin_2d() + diagramWidth_2d()))
      return false;
    if (y < topMargin_2d())
      return false;
    if (y > (topMargin_2d() + diagramHeight_2d()))
      return false;
    return true;
  }

  public int checkPopupClick(int x, int y)
  {
    float ix, z, n;

    if (in2dDiagram(x, y) == false)
      return 0;

    float t = getRetentionTime(x);
    float intensity = getIntensity(y);

    for (int i = 1; i < cr_.Value.length; i++) {
      if (cr_.Value[i][0] < t)
        continue;
      n = cr_.Value[i][0] - cr_.Value[i-1][0];
      if (raw_ == false) {
        z = cr_.Value[i][2] - cr_.Value[i-1][2];
        ix = cr_.Value[i-1][2] + z / n * (t - cr_.Value[i-1][0]);
      } else {
        z = cr_.Value[i][1] - cr_.Value[i-1][1];
        ix = cr_.Value[i-1][1] + z / n * (t - cr_.Value[i-1][0]);
      }
      if (intensity < ix)
        return 1;
      return 0;
    }
    return 0;
  }
  
  public boolean anyUserProbeThere(int x)
  {
    CgProbe cp = getUserProbe(x);
    if (cp!=null)
      return true;
    else
      return false;
      
  }
  
  public CgProbe getUserProbe(int x)
  {
    CgProbe cp;
    float t = getRetentionTime(x);
    for (int i = 0; i < this.selectedProbes_.size(); i++) {
      cp = (CgProbe)this.selectedProbes_.get(i);
      if (cp.LowerValley > t || cp.UpperValley < t)
        continue;
////      if (cp.Mz > (cr_.Mz + this.getThreshold())
////          || cp.Mz < (cr_.Mz - this.getThreshold()))
//      System.out.println("IsotopeNumbers: "+cp.isotopeNumber+";"+this.isotopeNumber_);
      if (cp.isotopeNumber!=this.isotopeNumber_)
        continue;
      return cp;
    }
    return null;
  }
  
  public boolean anySystemProbeThere(int x)
  {
    CgProbe cp = this.getSystemProbe(x);
    if (cp !=null)
      return true;
    else
      return false;
  }

  public CgProbe getSystemProbe(int x)
  {
    CgProbe cp;
    float t = getRetentionTime(x);
    for (int i = 0; i < this.storedProbes_.size(); i++) {
      cp = (CgProbe) ((this.storedProbes_).get(i));
      if (cp.LowerValley > t || cp.UpperValley < t)
        continue;
      if (cp.isotopeNumber!=this.isotopeNumber_)
        continue;
      boolean removed = false;
      for (int j = 0; j != this.removedStoredProbes_.size(); j++) {
        if (this.sameCgProbe((CgProbe) this.removedStoredProbes_.get(j), cp))
          removed = true;
      }
      if (!removed)
        return cp;
      else
        return null;
    }
    return null;
  }
  

  private void removeStoredArea(int x)
  {
    CgProbe cp;
    float t = getRetentionTime(x);
    for (int i = 0; i < this.storedProbes_.size(); i++) {
      cp = (CgProbe) ((storedProbes_).get(i));
      if (cp.LowerValley > t || cp.UpperValley < t)
        continue;
//      if (cp.Mz > (cr_.Mz + this.getThreshold())
//          || cp.Mz < (cr_.Mz - this.getThreshold()))
      if (cp.isotopeNumber!=this.isotopeNumber_)
        continue;
      this.removedStoredProbes_.add(cp);
      break;
    }
    this.draw2DDiagram(Float.valueOf(String.valueOf(currentStart_))/this.resolutionFactor_, Float.valueOf(String.valueOf(currentStop_))/this.resolutionFactor_, this.raw_);
    this.update(this.getGraphics());
    //    m_2dPnl.preparePaint();
//    Draw2DDiagram(m_2dPnl, m_raw);
//    this.refresh3dGraphics();
  }

  private boolean sameCgProbe(CgProbe probe1, CgProbe probe2)
  {
    boolean ok = false;
    /**** Here it is questionable if the CgProbe.Peak should be compared (for enhanced valley method) because the peak
     * can shift and then the area is painted green instead of red
     */
    if (probe1.Area >= (probe2.Area - 0.02)
        && probe1.Area <= (probe2.Area + 0.02)
        && probe1.AreaError >= (probe2.AreaError - 0.02)
        && probe1.AreaError <= (probe2.AreaError + 0.02)
        && probe1.Background >= (probe2.Background - 0.02)
        && probe1.Background <= (probe2.Background + 0.02)
        && probe1.Charge == probe2.Charge
        && probe1.LowerMzBand >= (probe2.LowerMzBand - this.getThreshold())
        && probe1.LowerMzBand <= (probe2.LowerMzBand + this.getThreshold())
        && probe1.LowerValley >= (probe2.LowerValley - 0.02)
        && probe1.LowerValley <= (probe2.LowerValley + 0.02)
        && probe1.Mz >= (probe2.Mz - this.getThreshold())
        && probe1.Mz <= (probe2.Mz + this.getThreshold())
        && probe1.Peak >= (probe2.Peak - 0.02)
        && probe1.Peak <= (probe2.Peak + 0.02)
        && probe1.UpperMzBand >= (probe2.UpperMzBand - this.getThreshold())
        && probe1.UpperMzBand <= (probe2.UpperMzBand + this.getThreshold())
        && probe1.UpperValley >= (probe2.UpperValley - 0.02)
        && probe1.UpperValley <= (probe2.UpperValley + 0.02)) {
      ok = true;
    }
    return ok;
  }
  
  private void deleteArea(int x)
  {
    float t = getRetentionTime(x);
    for (int i = 0; i < this.selectedProbes_.size(); i++) {
      CgProbe cp = (CgProbe) ((this.selectedProbes_).get(i));
      if (cp.LowerValley > t || cp.UpperValley < t)
        continue;
      if (cp.isotopeNumber!=this.isotopeNumber_)
//      if (cp.Mz > (cr_.Mz + this.getThreshold())
//          || cp.Mz < (cr_.Mz - this.getThreshold()))
        continue;
      this.selectedProbes_.remove(cp);
      break;
    }
    this.draw2DDiagram(Float.valueOf(String.valueOf(currentStart_))/this.resolutionFactor_, Float.valueOf(String.valueOf(currentStop_))/this.resolutionFactor_, this.raw_);
    this.update(this.getGraphics());
  }
  
  private void determineArea(int x, int method)
  {
    float t = getRetentionTime(x);
    CgProbe cx = null;
    if (method==LipidomicsDefines.Valley3DMethod){
      try {
        cx = analyzer_.detectPeakThreeD(cr_,LipidomicsAnalyzer.findIndexByTime(t, cr_),false,charge_,msLevel_);
      } catch (CgException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      } catch (QuantificationException qcx){
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(),"Warning", qcx.getMessage());
        return;
      }

    }else
      cx = LipidomicsAnalyzer.calculateOneArea(cr_, LipidomicsAnalyzer.findIndexByTime(t, cr_), method,charge_);
    cx.isotopeNumber = this.isotopeNumber_;
    if (cx.AreaStatus != CgAreaStatus.OK){
      @SuppressWarnings("unused")
      WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The selected peak is too small to be valid");
      return;
    }
    addProbeToTheCache(cx);
  }
  
  private void addProbeToTheCache(CgProbe cx){
    // First check if a removed system probe is there, then restore it
    if (!this.isRemovedSystemProbeThere(cx)){
      CgProbe overlappingProbe = overlappingToOtherProbe(cx);
      // Check if it overlaps
      if (overlappingProbe!=null){
        @SuppressWarnings("unused")
        WarningMessage dlg = new WarningMessage(new JFrame(), "Warning", "The selected probe area will overlap with the area "+overlappingProbe.Area+"; start: "+cx.LowerValley/60+"min stop: "+cx.UpperValley/60+"min");
      }else{
        this.selectedProbes_.add(cx);
      }
      
    }
    this.draw2DDiagram(Float.valueOf(String.valueOf(currentStart_))/this.resolutionFactor_, Float.valueOf(String.valueOf(currentStop_))/this.resolutionFactor_, this.raw_);
    this.update(this.getGraphics());
    
  }

  private boolean isRemovedSystemProbeThere(CgProbe newProbe){
    boolean isThere = false;
    int count = 0;
    int position = -1;
    for (CgProbe aProbe : this.removedStoredProbes_){
      if (this.sameCgProbe(newProbe,aProbe)){
        isThere = true;
        position = count;
      }
      count++;
    }
    if (isThere){
      this.removedStoredProbes_.remove(position);
    }
    return isThere;
  }
  
  private CgProbe overlappingToOtherProbe(CgProbe newProbe){
    CgProbe overlapProbe = null;
    for (CgProbe aProbe : this.storedProbes_){
      if ((newProbe.Mz > (aProbe.Mz - this.getThreshold()))
          && (newProbe.Mz < (aProbe.Mz + this.getThreshold()))) {
        if ((newProbe.UpperValley > aProbe.LowerValley && newProbe.LowerValley < aProbe.LowerValley)||
            (newProbe.LowerValley < aProbe.LowerValley && newProbe.UpperValley > aProbe.UpperValley)||
            (newProbe.LowerValley+1>aProbe.LowerValley && newProbe.UpperValley-1<aProbe.UpperValley) ||
            (newProbe.LowerValley-1<aProbe.LowerValley && newProbe.UpperValley+1>aProbe.UpperValley)){
          // this is just for testing
          boolean isRemovedOverlapProbe = false;
          for (CgProbe probe: this.removedStoredProbes_){
            if (probe.Area+1>aProbe.Area&&probe.Area-1<aProbe.Area)
              isRemovedOverlapProbe = true;
            
          }
//          System.out.println("Is removed overlap: "+isRemovedOverlapProbe);
          if (!isRemovedOverlapProbe)
            overlapProbe = aProbe;
        }
      }  
    }
    for (CgProbe aProbe : this.selectedProbes_){
      if ((newProbe.Mz > (aProbe.Mz - this.getThreshold()))
          && (newProbe.Mz < (aProbe.Mz + this.getThreshold()))) {
        if ((newProbe.UpperValley > aProbe.LowerValley && newProbe.LowerValley < aProbe.LowerValley)||
            (newProbe.LowerValley < aProbe.LowerValley && newProbe.UpperValley > aProbe.UpperValley)||
            (newProbe.LowerValley+1>aProbe.LowerValley && newProbe.UpperValley-1<aProbe.UpperValley) ||
            (newProbe.LowerValley-1<aProbe.LowerValley && newProbe.UpperValley+1>aProbe.UpperValley)){
          overlapProbe = aProbe;
        }
      }  
    }
    return overlapProbe;  
  }
  
  public void determineArea(int mode){
    this.determineArea(this.lastPopupX, mode);
  }
  
  public void deleteArea(){
    if (this.anySystemProbeThere(this.lastPopupX))
      removeStoredArea(this.lastPopupX);
    else
      deleteArea(this.lastPopupX);
  }
  
  public Vector<Vector<CgProbe>> getAllSelectedProbes(){
    Vector<Vector<CgProbe>> allProbes = new Vector<Vector<CgProbe>>();
    Vector<CgProbe> storedProbes = new Vector<CgProbe>();
    for (CgProbe probe : this.storedProbes_){
      boolean isRemoved = false;
      for (CgProbe removedProbe : this.removedStoredProbes_){
        if (this.sameCgProbe(probe, removedProbe)){
          isRemoved = true;
        }
      }
      if (!isRemoved){
        storedProbes.add(probe);
      }
    }
    allProbes.add(storedProbes);
    allProbes.add(this.selectedProbes_);
    return allProbes;
  }
  
  public void setStoredProbes(Vector<CgProbe> probes){
    this.storedProbes_ = probes;
    this.selectedProbes_ = new Vector<CgProbe>();
    this.removedStoredProbes_ = new Vector<CgProbe>();
    this.draw2DDiagram(Float.valueOf(String.valueOf(currentStart_))/this.resolutionFactor_, Float.valueOf(String.valueOf(currentStop_))/this.resolutionFactor_, raw_);
    this.update(this.getGraphics());
  }
  
  private CgProbe createManuallyAddedProbe(AreaSettingVO settings){
    CgProbe probeToAdd = null;
    double startTime = settings.getStartTime()*this.timeCorrectionFactor_;
    double stopTime = settings.getStopTime()*this.timeCorrectionFactor_;
//    boolean calculatePeakParametersCompletelyNew = false;
    if (probeToAdjust_!=null){
      // remove the old probe from the cache
      if (this.isProbeToAdjustSystemProbe_){
        removeStoredArea(this.lastPopupX);
      } else {
        this.deleteArea();
      }
//      if (probeToAdjust_ instanceof Probe3D){
//        Probe3D threeDProbe = (Probe3D)probeToAdjust_;
//        // The probe is a 3D probe and will remain one
//        if (settings.isThreeD()){
//        
//        // The 3D probe will become a 2D one
//        }else{
//          
//        }
//      }else{
//        // The probe is a 2D probe and should become a 3D probe 
//        if (settings.isThreeD()){
//        
//        // The probe is a 2D probe and will remain one
//        }else{
//          
//        }
//      }
      
    }//else{
//      calculatePeakParametersCompletelyNew = true;
//    }
//    if (calculatePeakParametersCompletelyNew){
    try{
    //TODO: the msLevel has to be put here somehow
      probeToAdd = analyzer_.getProbeByManualSettings(cr_, settings.isThreeD(), startTime, stopTime, (float)settings.getStartMz(), (float)settings.getStopMz(),charge_,msLevel_);

    }
    catch (CgException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }  
//    }
    return probeToAdd;
  }
  
  private static float getCorrespondingValue(float[] valueArray, float maxIntensity, boolean raw, boolean relativeValue){
    float value = 0;
    if (raw) value = valueArray[1];
    else value = valueArray[2];
    return getCorrespondingValue(value, maxIntensity, raw, relativeValue);
  }
  
  private static float getCorrespondingValue(float oldValue, float maxIntensity, boolean raw, boolean relativeValue){
    float value = oldValue;
    if (relativeValue) value = (value*100f)/maxIntensity;
    return value;
  }
  
  public void setRelativeValues(boolean relativeValues){
    this.relativeIntensity_ = relativeValues;
  }
  
  protected static int[] getCoordinatesInDiagram(float[][] chromValues, int pos, int x0, int y0, int w0, int h0,
      float m_minDispTime2d, float m_maxDispTime2d, float maxIntensity, float m_2dGain, boolean raw, boolean relativeValue,
      boolean stopAtH0){
    int[] coords = new int[2];
    coords[0] = x0 + (int) (w0 * (chromValues[pos][0] - m_minDispTime2d) / (m_maxDispTime2d - m_minDispTime2d));
    float maxInt = getMaxInt(maxIntensity,relativeValue);
    float value = getCorrespondingValue(chromValues[pos],maxIntensity,raw,relativeValue);
    coords[1] = y0 - (int) ((h0 * value) * m_2dGain / maxInt);
    if (stopAtH0 && coords[1]<(y0-h0)) coords[1] = (y0-h0);
    return coords;
  }
  
  private static float getMaxInt(float maxIntensity, boolean relativeValue){
    float maxInt = 100f;
    if (!relativeValue) maxInt = maxIntensity;
    maxInt = maxInt*defaultIntensityRangeMultiplicator_;    
    return maxInt;
  }
  
  public void paintMs2Position(Float ms2Position){
    this.ms2Position_ = ms2Position;
//    draw2DDiagram(chromMzStart_, chromMzStop_, raw_);
    repaint();
  }
  
  
  public float getMinDispTime2D()
  {
    return m_minTime_;
  }
  
  public float getMaxDispTime2D()
  {
    return m_maxTime_;
  }
  
  public float getM2dGain()
  {
    return m_2dGain_;
  }
  
  public void setM2dGain(float m_2dGain)
  {
    this.m_2dGain_ = m_2dGain;
  }
  
  public void nextSpectrum(){}
  public void nextSpectrum(LipidParameterSet param, Vector<RangeColor> rangeColors){}
  public void previousSpectrum(){}
  public void previousSpectrum(LipidParameterSet param, Vector<RangeColor> rangeColors){}
  public void refresh(LipidParameterSet param, Vector<RangeColor> rangeColors){}
  public String getSpectSelectedText(){return null;}
  public String getRtSelectedText(){return null;}
  public String getPrecursorMassSelected(){return null;}
  public float[] getRTRange(){return new float[2];}
  public int getSpectrumSelected(){return -1;};
  public void clearRangeColors(){};
  public void setAnnotationThreshold(double cutoff){};
  
  public BufferedImage getImage(){
    BufferedImage chart=new BufferedImage(getWidth(),getHeight(),
        BufferedImage.TYPE_3BYTE_BGR);
    draw2DDiagram(chart.getGraphics());
    return chart;
  }

}
