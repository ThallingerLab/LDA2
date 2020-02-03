package at.tugraz.genome.vos;

import java.util.LinkedHashMap;

public class Mix1Standards
{
  
  public static LinkedHashMap<String,ReferenceInfoVO> getSphBaseStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d17:0", new ReferenceInfoVO("d17:0",new double[]{2.6d},true,false));
    standards.get("d17:0").addAdduct("H", 288.289705857d);
    standards.get("d17:0").addAdduct("-OH", 270.279141157d);
    standards.get("d17:0").addAdduct("Na", 310.2716515577d);
    standards.put("d17:1", new ReferenceInfoVO("d17:1",new double[]{2.3d},true,false));
    standards.get("d17:1").addAdduct("H", 286.274055787d);
    standards.get("d17:1").addAdduct("-OH", 268.263491087d);
    standards.get("d17:1").addAdduct("Na", 308.256001487d);
    standards.put("t18:0", new ReferenceInfoVO("t18:0",new double[]{2.6d},true,false));
    standards.get("t18:0").addAdduct("H", 318.300270557d);
    standards.get("t18:0").addAdduct("-OH", 300.289705857d);
    standards.get("t18:0").addAdduct("Na", 340.282216257d);
    standards.put("m18:1", new ReferenceInfoVO("m18:1",new double[]{3.9d},false,false));
    standards.get("m18:1").addAdduct("H", 318.300270557d);
    standards.get("m18:1").addAdduct("-OH", 300.289705857d);
    standards.get("m18:1").addAdduct("Na", 340.282216257d);
    
    return standards;
  }
  
  public static LinkedHashMap<String,ReferenceInfoVO> getCerStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d16:1", new ReferenceInfoVO("d14:1/n2:0",new double[]{1.4d},true,false));
    standards.get("d16:1").addAdduct("-H", 284.223117677d);
    standards.get("d16:1").addAdduct("HCOO", 330.228597007d);
    standards.get("d16:1").addAdduct("Cl", 320.199795977d);
    standards.get("d16:1").addAdduct("H", 286.237670277d);
    standards.get("d16:1").addAdduct("-OH", 268.227105577d);
    standards.get("d16:1").addAdduct("Na", 308.219615977d);
    standards.put("t26:0", new ReferenceInfoVO("t18:0/n8:0",new double[]{9.9d},true,false));
    standards.get("t26:0").addAdduct("-H", 442.390183077d);
    standards.get("t26:0").addAdduct("HCOO", 488.395662407d);
    standards.get("t26:0").addAdduct("Cl", 478.366861377d);
    standards.get("t26:0").addAdduct("H", 444.404735677d);
    standards.get("t26:0").addAdduct("-OH", 426.394170977d);
    standards.get("t26:0").addAdduct("Na", 466.386681377d);
    standards.put("m30:0", new ReferenceInfoVO("m18:0/n12:0",new double[]{19.8d},false,false));
    standards.get("m30:0").addAdduct("-H", 466.462954097d);
    standards.get("m30:0").addAdduct("HCOO", 512.468433427d);
    standards.get("m30:0").addAdduct("Cl", 502.439632397d);
    standards.get("m30:0").addAdduct("H", 468.477506697d);
    standards.get("m30:0").addAdduct("-OH", 450.466941997d);
    standards.get("m30:0").addAdduct("Na", 490.459452397d);
    standards.put("d30:1", new ReferenceInfoVO("d18:1/n12:0",new double[]{17.2d},true,false));
    standards.get("d30:1").addAdduct("-H", 480.442218657d);
    standards.get("d30:1").addAdduct("HCOO", 526.447697987d);
    standards.get("d30:1").addAdduct("Cl", 516.418896957d);
    standards.get("d30:1").addAdduct("H", 482.456771257d);
    standards.get("d30:1").addAdduct("-OH", 464.446206557d);
    standards.get("d30:1").addAdduct("Na", 504.438716957d);
    standards.put("t30:1", new ReferenceInfoVO("d18:1/h12:0",new double[]{15.9d},true,false));
    standards.get("t30:1").addAdduct("-H", 496.437133287d);
    standards.get("t30:1").addAdduct("HCOO", 542.442612617d);
    standards.get("t30:1").addAdduct("Cl", 532.413811587d);
    standards.get("t30:1").addAdduct("H", 498.451685887d);
    standards.get("t30:1").addAdduct("-OH", 480.441121187d);
    standards.get("t30:1").addAdduct("Na", 520.433631587d);
    standards.put("d32:1", new ReferenceInfoVO("d16:1/n16:0",new double[]{19.7d},true,false));
    standards.get("d32:1").addAdduct("-H", 508.473518797d);
    standards.get("d32:1").addAdduct("HCOO", 554.478998127d);
    standards.get("d32:1").addAdduct("Cl", 544.450197097d);
    standards.get("d32:1").addAdduct("H", 510.488071397d);
    standards.get("d32:1").addAdduct("-OH", 492.477506697d);
    standards.get("d32:1").addAdduct("Na", 532.470017097d);  
    standards.put("t34:0", new ReferenceInfoVO("t18:0/n16:0",new double[]{20.8d},true,false));
    standards.get("t34:0").addAdduct("-H", 554.515383637d);
    standards.get("t34:0").addAdduct("HCOO", 600.520862967d);
    standards.get("t34:0").addAdduct("Cl", 590.492061937d);
    standards.get("t34:0").addAdduct("H", 556.529936237d);
    standards.get("t34:0").addAdduct("-OH", 538.519371537d);
    standards.get("t34:0").addAdduct("Na", 578.511881937d);    
    standards.put("d34:2", new ReferenceInfoVO("d18:2/n16:0",new double[]{21.0d},true,false));
    standards.get("d34:2").addAdduct("-H", 534.489168867d);
    standards.get("d34:2").addAdduct("HCOO", 580.494648197d);
    standards.get("d34:2").addAdduct("Cl", 570.465847167d);
    standards.get("d34:2").addAdduct("H", 536.503721467d);
    standards.get("d34:2").addAdduct("-OH", 518.493156767d);
    standards.get("d34:2").addAdduct("Na", 558.485667167d);    
    standards.put("d36:0", new ReferenceInfoVO("d18:0/n18:0",new double[]{24.9d},true,false));
    standards.get("d36:0").addAdduct("-H", 566.551769147d);
    standards.get("d36:0").addAdduct("HCOO", 612.557248477d);
    standards.get("d36:0").addAdduct("Cl", 602.528447447d);
    standards.get("d36:0").addAdduct("H", 568.566321747d);
    standards.get("d36:0").addAdduct("-OH", 550.555757047d);
    standards.get("d36:0").addAdduct("Na", 590.548267447d);
    standards.put("t36:0", new ReferenceInfoVO("d18:0/h18:0",new double[]{23.4d,23.9d},true,false));
    standards.get("t36:0").addAdduct("-H", 582.546683777d);
    standards.get("t36:0").addAdduct("HCOO", 628.552163107d);
    standards.get("t36:0").addAdduct("Cl", 618.523362077d);
    standards.get("t36:0").addAdduct("H", 584.561236377d);
    standards.get("t36:0").addAdduct("-OH", 566.550671677d);
    standards.get("t36:0").addAdduct("Na", 606.543182077d);    
    standards.put("d36:2", new ReferenceInfoVO("d18:1/n18:1",new double[]{22.7d},true,false));
    standards.get("d36:2").addAdduct("-H", 562.520469007d);
    standards.get("d36:2").addAdduct("HCOO", 608.525948337d);
    standards.get("d36:2").addAdduct("Cl", 598.497147307d);
    standards.get("d36:2").addAdduct("H", 564.535021607d);
    standards.get("d36:2").addAdduct("-OH", 546.524456907d);
    standards.get("d36:2").addAdduct("Na", 586.516967307d);
    standards.put("t42:0", new ReferenceInfoVO("t18:0/n24:0",new double[]{28.8d},true,false));
    standards.get("t42:0").addAdduct("-H", 666.640584197d);
    standards.get("t42:0").addAdduct("HCOO", 712.646063527d);
    standards.get("t42:0").addAdduct("Cl", 702.617262497d);
    standards.get("t42:0").addAdduct("H", 668.655136797d);
    standards.get("t42:0").addAdduct("-OH", 650.644572097d);
    standards.get("t42:0").addAdduct("Na", 690.637082497d);    
    standards.put("q42:0", new ReferenceInfoVO("t18:0/h24:0",new double[]{28.1d},true,false));
    standards.get("q42:0").addAdduct("-H", 682.635498827d);
    standards.get("q42:0").addAdduct("HCOO", 728.640978157d);
    standards.get("q42:0").addAdduct("Cl", 718.612177127d);
    standards.get("q42:0").addAdduct("H", 684.650051427d);
    standards.get("q42:0").addAdduct("-OH", 666.639486727d);
    standards.get("q42:0").addAdduct("Na", 706.631997127d);    
    standards.put("m42:1", new ReferenceInfoVO("m18:0/n24:1",new double[]{29.6d},false,false));
    standards.get("m42:1").addAdduct("-H", 632.635104867d);
    standards.get("m42:1").addAdduct("HCOO", 678.640584197d);
    standards.get("m42:1").addAdduct("Cl", 668.611783167d);
    standards.get("m42:1").addAdduct("H", 634.649657467d);
    standards.get("m42:1").addAdduct("-OH", 616.639092767d);
    standards.get("m42:1").addAdduct("Na", 656.631603167d);    
    standards.put("d43:1", new ReferenceInfoVO("d18:1/n25:0",new double[]{30.4d},true,false));
    standards.get("d43:1").addAdduct("-H", 662.645669567);
    standards.get("d43:1").addAdduct("HCOO", 708.651148897d);    
    standards.get("d43:1").addAdduct("Cl", 698.622347867d);
    standards.get("d43:1").addAdduct("H", 664.660222167d);
    standards.get("d43:1").addAdduct("-OH", 646.649657467d);
    standards.get("d43:1").addAdduct("Na", 686.6421678677d);
    return standards;
  }
  
  public static LinkedHashMap<String,ReferenceInfoVO> getCer1PStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d30:1", new ReferenceInfoVO("d18:1/n12:0",new double[]{13.2d},true,false));
    standards.get("d30:1").addAdduct("-H", 560.408549582);
    standards.get("d30:1").addAdduct("H", 562.423102182d);
    standards.get("d30:1").addAdduct("-OH", 544.412537482d);
    standards.get("d30:1").addAdduct("Na", 584.405047882d);
    standards.put("d34:0", new ReferenceInfoVO("d18:0/n16:0",new double[]{20.4d},true,false));
    standards.get("d34:0").addAdduct("-H", 618.486799932);
    standards.get("d34:0").addAdduct("H", 620.501352532d);
    standards.get("d34:0").addAdduct("-OH", 602.490787832d);
    standards.get("d34:0").addAdduct("Na", 642.483298232d);
    return standards;
  }
    
  
  public static LinkedHashMap<String,ReferenceInfoVO> getHexCerStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d26:1", new ReferenceInfoVO("d18:1/n8:0",new double[]{9.4d},true,false));
    standards.get("d26:1").addAdduct("-H", 586.432441877d);
    standards.get("d26:1").addAdduct("HCOO", 632.437921207d);
    standards.get("d26:1").addAdduct("Cl", 622.409120177d);
    standards.get("d26:1").addAdduct("H", 588.446994477d);
    standards.get("d26:1").addAdduct("-OH", 570.436429777d);
    standards.get("d26:1").addAdduct("Na", 610.428940177d);
    standards.put("d30:1", new ReferenceInfoVO("d18:1/n12:0",new double[]{14.9d},true,false));
    standards.get("d30:1").addAdduct("-H", 642.495042157d);
    standards.get("d30:1").addAdduct("HCOO", 688.500521487d);
    standards.get("d30:1").addAdduct("Cl", 678.471720457d);
    standards.get("d30:1").addAdduct("H", 644.509594757d);
    standards.get("d30:1").addAdduct("-OH", 626.499030057d);
    standards.get("d30:1").addAdduct("Na", 666.491540457d);
    standards.put("d34:0", new ReferenceInfoVO("d18:0/n16:0",new double[]{20.7d},true,false));
    standards.get("d34:0").addAdduct("-H", 700.573292507d);
    standards.get("d34:0").addAdduct("HCOO", 746.578771837d);
    standards.get("d34:0").addAdduct("Cl", 736.549970807d);
    standards.get("d34:0").addAdduct("H", 702.587845107d);
    standards.get("d34:0").addAdduct("-OH", 684.577280407d);
    standards.get("d34:0").addAdduct("Na", 724.569790807d);
    standards.put("d36:1", new ReferenceInfoVO("d18:1/n18:0",new double[]{22.3d},true,false));
    standards.get("d36:1").addAdduct("-H", 726.588942577d);
    standards.get("d36:1").addAdduct("HCOO", 772.594421907d);
    standards.get("d36:1").addAdduct("Cl", 762.565620877d);
    standards.get("d36:1").addAdduct("H", 728.603495177d);
    standards.get("d36:1").addAdduct("-OH", 710.592930477d);
    standards.get("d36:1").addAdduct("Na", 750.585440877d);
    standards.put("t36:1", new ReferenceInfoVO("d18:1/h18:0",new double[]{21.5d},true,false));
    standards.get("t36:1").addAdduct("-H", 742.583857207d);
    standards.get("t36:1").addAdduct("HCOO", 788.589336537d);
    standards.get("t36:1").addAdduct("Cl", 778.560535507d);
    standards.get("t36:1").addAdduct("H", 744.598409807d);
    standards.get("t36:1").addAdduct("-OH", 726.587845107d);
    standards.get("t36:1").addAdduct("Na", 750.585440877d);
    return standards;
  }
  
  
  public static LinkedHashMap<String,ReferenceInfoVO> getSMStandards(){
    LinkedHashMap<String,ReferenceInfoVO> standards = new LinkedHashMap<String,ReferenceInfoVO>();
    // the first boolean is if the standard shall be used in the evaluation, and the second if a position is detectable
    standards.put("d30:1", new ReferenceInfoVO("d18:1/n12:0",new double[]{13.5d},true,false));
    standards.get("d30:1").addAdduct("HCOO", 691.503178299d);
    standards.get("d30:1").addAdduct("Cl", 681.474377269d);
    standards.get("d30:1").addAdduct("H", 647.512251569d);
    standards.get("d30:1").addAdduct("Na", 669.494197269d);
    standards.put("d36:2", new ReferenceInfoVO("d18:1/n18:1",new double[]{19.6d},true,false));
    standards.get("d36:2").addAdduct("HCOO", 773.581428649d);
    standards.get("d36:2").addAdduct("Cl", 763.552627619d);
    standards.get("d36:2").addAdduct("H", 729.590501919d);
    standards.get("d36:2").addAdduct("Na", 751.572447619d);
    standards.put("d42:2", new ReferenceInfoVO("d18:1/n24:1",new double[]{25.7d},true,false));
    standards.get("d42:2").addAdduct("HCOO", 857.675329069d);
    standards.get("d42:2").addAdduct("Cl", 847.646528039d);
    standards.get("d42:2").addAdduct("H", 813.684402339d);
    standards.get("d42:2").addAdduct("Na", 835.666348039d);
    return standards;
  }
  
}
