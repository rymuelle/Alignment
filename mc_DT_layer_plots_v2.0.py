import ROOT as r
import struct, math, os, sys
from array import array

selectedDT = (0,4,10) # wheel station chamber  
specialSuffix = ""

# suppress "Info in..." messages, as there will be a lot
r.gErrorIgnoreLevel = r.kWarning


nSigma = 1.5 # number of sigma to do gaussian fit with
ME13pinR = 595.1500244141 # radius of ME13 ring alignment pin
DT_0_4_xpinR = 700 # approx radius of DT sector 4 alignment pin
isMC = False
count_max = 250000*60
count_max = -1

fid_cut_x_1 = 170
fid_cut_x_2 = -170
#fid_cut_x_1 = 100
#fid_cut_x_2 = -100
fid_cut_y_1 = 170
fid_cut_y_2 = -170
#fid_cut_y_1 = 100
#fid_cut_y_2 = -100

print ">>> Args:", sys.argv
if(len(sys.argv) > 1):
    filename = sys.argv[1]
    dtStr = sys.argv[1].split("DT")[1]
    dtStr = dtStr.split(".root")[0]
    dtInfo = dtStr.split("_")
    #selectedDT = ( dtInfo[0], int(dtInfo[1]), int(dtInfo[2]), int(dtInfo[3]) )
#    selectedDT = (  int(dtInfo[1]), int(dtInfo[2]), int(dtInfo[3]) )
    selectedDT = ( 0, 4, 1)
    if(dtInfo[4] == "MC"):
        specialSuffix = dtInfo[5]
        isMC = True
    else:
        specialSuffix = dtInfo[4]
    isMC = True

    print ">>> Found argument with filename:", filename
    print ">>> Setting selectedDT to", selectedDT
else:
    print ">>> Using default file:", filename
    print ">>> Setting selectedDT to default:", selectedDT

if(specialSuffix is not ""): print ">>> Special suffix is",specialSuffix
    
prefix = "DT_%d_%d_%d_n_%i_track_cuts" % ( selectedDT[0], selectedDT[1], selectedDT[2], count_max)
if(isMC): prefix += "_MC"
prefix += "_%s_1.5_sigma_fit_plots/" % (specialSuffix)
print ">>> Output into folder",prefix


def rProjection(localr, angle):
    return (localr + ME13pinR) / math.cos(angle) - ME13pinR


def fitCut(hist, sigmas, opts):
    lower, upper = hist.GetMean()-sigmas*hist.GetRMS(), hist.GetMean()+sigmas*hist.GetRMS()
    hist.Fit("gaus",opts, "", lower, upper)

def getFitParams(hist):
    fit = hist.GetFunction("pol1")
    
    p0 = fit.GetParameter(0) # offset from x axis
    p0e = fit.GetParError(0) # offset from x axis
    
    p1 = fit.GetParameter(1) # slope
    p1e = fit.GetParError(1) # slope
    
    return p0, p0e, p1, p1e

def evalInCenter(hist):
    fit = hist.GetFunction("pol1")
    p0 = fit.GetParameter(0) # offset from x axis
    p0e = fit.GetParError(0) # offset from x axis
    p1 = fit.GetParameter(1) # slope
    p1e = fit.GetParError(1) # slope

    # evaluate at layer 3.5
    val = p0+3.5*p1
    err = math.sqrt(p0e**2 + 3.5**2 * p1e**2)
    
    return val, err


def getFitParamsGauss(hist):
    fit = hist.GetFunction("gaus")
    
    p0 = fit.GetParameter(0) # const
    p0e = fit.GetParError(0) # const
    
    p1 = fit.GetParameter(1) # mean
    p1e = fit.GetParError(1) # mean error
    
    #print p0, p0e, p1, p1e
    return p0, p0e, p1, p1e

def fixFlooredBins(hist,minZ=-3.5):
    for x in range(hist.GetNbinsX()+1):
        for y in range(hist.GetNbinsY()+1):
            mu = hist.GetBinContent(x,y)
            if(mu <= minZ):
                hist.SetBinContent(x,y,minZ)

def zoomXrange(h):
    nonEmptyBins = [(e,h.GetXaxis().GetBinCenter(e)) for e in range(h.GetNbinsX()) if h.GetBinContent(e)>0] 
    minX, maxX = min(nonEmptyBins)[1], max(nonEmptyBins)[1]
    tolerance = (maxX-minX)*0.15
    h.GetXaxis().SetRangeUser(minX-tolerance,maxX+tolerance)

def binData(indexPt, Nbins, minPt):
    d = {}

    binWidth = 1.0/abs(minPt)/Nbins
    
    for i in range(Nbins+1):
        d[i] = []
        
    for elem in indexPt:
        idx = elem[0]
        pt = elem[1]
        binIdx = int(math.floor(1./abs(pt)/binWidth));
        #print binIdx, pt
        d[binIdx].append(idx)
    
    return d
    
def getGoodMuonIndices(posIndexPt, negIndexPt, Nbins, minPt):
    dP = binData(posIndexPt, Nbins, minPt)
    dN = binData(negIndexPt, Nbins, minPt)
    indicesToConsider = []
#    for bin in range(Nbins):
#        minMuonsInBin = min(len(dP[bin]), len(dN[bin]))
#        indicesToConsider.extend(dP[bin][:minMuonsInBin])
#        indicesToConsider.extend(dN[bin][:minMuonsInBin])
#        
#        #print bin, len(dP[bin]), len(dN[bin]), minMuonsInBin #, d[bin]
#    return indicesToConsider

def applyRotation(x,y,phiz):
    xp = math.cos(phiz)*x - math.sin(phiz)*y
    yp = math.sin(phiz)*x + math.cos(phiz)*y
    return xp,yp

c1 = r.TCanvas("Canvas1", "Alignment Visualizations")

typePrefix = "DT %d/%d/%d" % ( selectedDT[0], selectedDT[1], selectedDT[2])
if(isMC):
    typePrefix += " #font[2]{#color[4]{MC}} "
else:
    typePrefix += " #font[2]{#color[2]{DATA}} "

nXbins, nYbins = 50, 50
XMax, YMax = 225., 225.
xyRatio = YMax/XMax 
#
resRange = 10
NUMLAYERS_X = 8
NUMLAYERS_Y = 8
#
#h2D_cnt_actual = []
#h2D_cnt_tracks = []
h2D_res_x_actual = []
h2D_res_x_tracks = []
h2D_res_y_actual = []
h2D_res_y_tracks = []
#h2D_pull_tracks = []
#h2D_pull_tracks = []
#h1D_res_x_rproj = []
##h1D_res_x_xproj = []
#h1D_res_x_rphiproj = []
#h1D_res_rphi_yproj = []
h1D_res_x_actual_x = []
h1D_res_x_track_x = []
h1D_res_x_track_x_2nd = []
h1D_res_x_track_x_1st = []
h1D_res_x_track_y = []
h1D_pull_actual_x = []
h1D_pull_tracks_x = []
h1D_res_x_actual_x_TH2F = []
h1D_res_x_track_x_TH2F = []
h1D_res_x = []
h1D_res_x_center_cut = []
h1D_res_x_top_cut = []
h1D_res_x_bottom_cut = []
h1D_res_y = []

n_x_chamber_cut = 20
h1D_res_x_cuts = []
h1D_res_x_cuts_sym = []
h1D_adjusted_res_x_cuts = []
h1D_adjusted_res_x_cuts_sym = []
tProfile_cut_bins = 100
for i in range(n_x_chamber_cut):
    title = "x Residual from %d to %d;x residual (cm);counts" % (i*2*XMax/n_x_chamber_cut-XMax, (i+1)*2*XMax/n_x_chamber_cut-XMax)
    h1D_res_x_cuts.append( r.TH1F("h1D_res_x_cuts_" + str(i), typePrefix + title, 200, -resRange, resRange) )
    h1D_res_x_cuts_sym.append( r.TH1F("h1D_res_x_cuts_sym_" + str(i), typePrefix + title, 200, -resRange, resRange) )
    title = "adjusted x Residual from %d to %d;adjusted x residual (cm);counts" % (i*2*XMax/n_x_chamber_cut-XMax, (i+1)*2*XMax/n_x_chamber_cut-XMax)
    h1D_adjusted_res_x_cuts.append( r.TH1F("h1D_adjusted_res_x_cuts_" + str(i), typePrefix + title, 200, -resRange, resRange) )
    h1D_adjusted_res_x_cuts_sym.append( r.TH1F("h1D_adjusted_res_x_cuts_sym_" + str(i), typePrefix + title, 200, -resRange, resRange) )

h1D_res_x_cuts.append( r.TH1F("h1D_res_x_cuts_" + str(i+1), typePrefix +"x Residual layer 1 ;x residual (cm);counts", 100, -resRange, resRange) )
h1D_adjusted_res_x_cuts.append( r.TH1F("h1D_adjusted_res_x_cuts_" + str(i+1), typePrefix +"adjusted x Residual layer 1 ;x residual (cm);counts", 100, -resRange, resRange) )
h1D_res_x_cuts_sym.append( r.TH1F("h1D_res_x_cuts_sym_" + str(i+1), typePrefix +"x Residual layer 1 ;x residual (cm);counts", 100, -resRange, resRange) )
h1D_adjusted_res_x_cuts_sym.append( r.TH1F("h1D_adjusted_res_x_cuts_sym_" + str(i+1), typePrefix +"adjusted x Residual layer 1 ;x residual (cm);counts", 100, -resRange, resRange) )

#h2D_res_x_mean_v_hit_x = r.TH2F("h2D_res_x_mean_v_hit_x", "x residual mean vs. track x; local x (cm); x residual (cm)", tProfile_cut_bins, -XMax, XMax,100, -resRange, resRange)
h2D_res_x_mean_v_hit_x = r.TH1F("h2D_res_x_mean_v_hit_x",  typePrefix + "x residual mean vs. track x; local x (cm); x residual mean (cm)", tProfile_cut_bins, -XMax, XMax)
h2D_res_x_median_v_hit_x = r.TH1F("h2D_res_x_median_v_hit_x",  typePrefix + "x residual median vs. track x; local x (cm); x residual median (cm)", tProfile_cut_bins, -XMax, XMax)
h2D_adjusted_res_x_mean_v_hit_x = r.TH1F("h2D_adjusted_res_x_mean_v_hit_x",  typePrefix + "adjusted x residual mean vs. track x; local x (cm);adjusted x residual mean (cm)", tProfile_cut_bins, -XMax, XMax)
h2D_adjusted_res_x_median_v_hit_x = r.TH1F("h2D_adjusted_res_x_median_v_hit_x",  typePrefix + "adjusted x residual median vs. track x; local x (cm);adjusted x residual median (cm)", tProfile_cut_bins, -XMax, XMax)
#h2D_res_x_RMS_v_hit_x = r.TProfile("h2D_res_x_RMS_v_hit_x", "x residual RMS vs. track x; local x (cm); x residual (cm)", tProfile_cut_bins, -XMax, XMax, -resRange, resRange)
h2D_res_x_RMS_v_hit_x = r.TH1F("h2D_res_x_RMS_v_hit_x", typePrefix + "x residual RMS vs. track x; local x (cm); x residual RMS (cm)", tProfile_cut_bins, -XMax, XMax)
h2D_res_x_skewness_v_hit_x = r.TH1F("h2D_res_x_skewness_v_hit_x", typePrefix + "x residual skewness vs. track x; local x (cm); skewness", tProfile_cut_bins, -XMax, XMax)
h2D_res_x_fit_p0_v_hit_x = r.TH1F("h2D_res_x_fit_p0_v_hit_x",  typePrefix + "x residual fit_p0 vs. track x; local x (cm); p0 (mm?)", tProfile_cut_bins, -XMax, XMax)
h2D_res_x_fit_p1_v_hit_x = r.TH1F("h2D_res_x_fit_p1_v_hit_x",  typePrefix + "x residual fit_p1 vs. track x; local x (cm); p1 (cm?)", tProfile_cut_bins, -XMax, XMax)

h2D_adjusted_res_x_RMS_v_hit_x = r.TH1F("h2D_adjusted_res_x_RMS_v_hit_x", typePrefix + "x  adjusted residual RMS vs. track x; local x (cm); x  adjusted residual RMS (cm)", tProfile_cut_bins, -XMax, XMax)
h2D_adjusted_res_x_skewness_v_hit_x = r.TH1F("h2D_adjusted_res_x_skewness_v_hit_x", typePrefix + "x  adjusted residual skewness vs. track x; local x (cm); skewness", tProfile_cut_bins, -XMax, XMax)
h2D_adjusted_res_x_fit_p0_v_hit_x = r.TH1F("h2D_adjusted_res_x_fit_p0_v_hit_x",  typePrefix + "x  adjusted residual fit_p0 vs. track x; local x (cm); p0 (mm?)", tProfile_cut_bins, -XMax, XMax)
h2D_adjusted_res_x_fit_p1_v_hit_x = r.TH1F("h2D_adjusted_res_x_fit_p1_v_hit_x",  typePrefix + "x  adjusted residual fit_p1 vs. track x; local x (cm); p1 (cm?)", tProfile_cut_bins, -XMax, XMax)
#h2D_statsig_tracks = []
#
##h1D_res_x_avg = r.TH1F("h1D_res_x_avg", "avg x residuals", 100,-8.0,8.0) 
#h1D_actual_localy = r.TH1F("h1D_actual_localy", typePrefix+"distribution of actual y hits (on L3);cm;counts",  48*2,-90,90)
#h1D_tracks_localy = r.TH1F("h1D_tracks_localy", typePrefix+"distribution of track y positions (on L3);cm;counts",  nYbins,-90,90)
#h1D_actual_angle = r.TH1F("h1D_actual_angle", typePrefix+"angular distribution of hits (on L3);phi;counts",  nYbins,-0.08,0.08)
#h1D_tracks_angle = r.TH1F("h1D_tracks_angle", typePrefix+"angular distribution of tracks (on L3);phi;counts",  nYbins,-0.08,0.08)
#h1D_rot_dxdr_layers = r.TH1F("h1D_rot_dxdr_layers", typePrefix+"phiz (dx/dr) vs layer;layer;dx/dr (urad)",  6,0.5,6.5  )
#h1D_trans_layers = r.TH1F("h1D_trans_layers", typePrefix+"x offset vs layer;layer; x offset (microns)",  6,0.5,6.5  )
#h2D_nlayers_hit = r.TProfile2D("h2D_nlayers_hit", typePrefix+"num layers hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
#h2D_nDT_hit = r.TProfile2D("h2D_nDT_hit", typePrefix+"num DTs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
#h2D_nDT_hit = r.TProfile2D("h2D_nDT_hit", typePrefix+"num DTs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
#h2D_nDTDT_hit = r.TProfile2D("h2D_nDTDT_hit", typePrefix+"num DTs+DTs hit;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
#h2D_nTracker_hit = r.TProfile2D("h2D_nTracker_hit", typePrefix+"num tracker hits;actual hit x;actual hit y",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)
h1D_pt = r.TH1F("h1D_pt", typePrefix+"p_{T};p_{T} (GeV/c);counts",  100, 20, 210)
h1D_pz = r.TH1F("h1D_pz", typePrefix+"p_{z};p_{z} (GeV/c);counts",  100, 20, 210)
h1D_p = r.TH1F("h1D_p", typePrefix+"p;p (GeV/c);counts",  100, 20, 210) 
h1D_eta = r.TH1F("h1D_eta", typePrefix+"#eta;#eta;counts", 800, -2.2, 2.2) 
h1D_phi = r.TH1F("h1D_phi", typePrefix+"#phi;#phi;counts",  800, -math.pi, math.pi)
h1D_x_occupancy = r.TH1F("h1D_x_occupancy", typePrefix+"h1D_x_occupancy;local x (cm) counts", 100, -XMax,XMax)
h1D_x_occupancy_whole_ring_special_chambers = r.TH1F("h1D_x_occupancy_whole_ring_special_chambers", typePrefix+"h1D_x_occupancy_whole_ring_special_chambers;local x (cm); counts", 1000, -100,6000)
h1D_x_occupancy_whole_ring = r.TH1F("h1D_x_occupancy_whole_ring", typePrefix+"h1D_x_occupancy_whole_ring;local x (cm) ;counts", 1000, -100,6000)

h2D_hit_x_res_x_whole_ring_special_chambers = r.TH2F("h2D_hit_x_res_x_whole_ring_special_chambers", typePrefix+"h2D_hit_x_res_x_whole_ring;local x (cm);x residual (cm)", 1000, -100, 6000, 100, -10, 10)
h2D_hit_x_res_x_whole_ring = r.TH2F("h2D_hit_x_res_x_whole_ring", typePrefix+"h2D_hit_x_res_x_whole_ring;local x (cm);x residual (cm)", 1000, -100,6000, 100, -resRange,resRange)

h2D_res_x_phi_eta = r.TH2F("h2D_res_x_phi_eta", typePrefix+"avg x res at actual hit positions;#phi;#eta",  500, -3.3,3.3,  500, -1.1,1.1)  
h2D_res_y_phi_eta = r.TH2F("h2D_res_y_phi_eta", typePrefix+"avg y res at actual hit positions;#phi;#eta",  500, -3.3,3.3,  500, -1.1,1.1)  
h2D_resx_res_y = r.TH2F("h2D_resx_res_y", typePrefix+"x res vs. y residual in their first layers;x res;y res", 100, -20,20, 100,-20,20)
#
#
#for i in range(NUMLAYERS):
for i in range(NUMLAYERS_X):
    laypfx = " #font[2]{L" + str(i+1) + "} "
    print ">>> Booking histos for layer",i+1
#    h2D_cnt_actual.append(  r.TH2F(laypfx + "h2D_cnt_actual", typePrefix+laypfx + "actual hit locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
#    h2D_cnt_tracks.append(  r.TH2F(laypfx + "h2D_cnt_tracks", typePrefix+laypfx + "track locations;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
#
    h2D_res_x_actual.append(  r.TProfile2D(laypfx + "h2D_res_x_actual", typePrefix+laypfx + "avg x res at actual hit positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
    h2D_res_x_tracks.append(  r.TProfile2D(laypfx + "h2D_res_x_tracks", typePrefix+laypfx + "avg x res at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
#
#    h2D_pull_tracks.append(  r.TProfile2D(laypfx + "h2D_pull_tracks", typePrefix+laypfx + "pull at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
#    h2D_statsig_tracks.append(  r.TH2F(laypfx + "h2D_statsig_tracks", typePrefix+laypfx + "significance at track positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
#
    h1D_res_x.append(  r.TH1F(laypfx + "h1D_res_x", typePrefix+laypfx + "x residuals;x residual (cm);counts",  100,-resRange,resRange)  )
    h1D_res_x_center_cut.append(  r.TH1F(laypfx + "h1D_res_x_center_cut", typePrefix+laypfx + "x residuals for center hit x;x residual (cm);counts",  100,-resRange,resRange)  )
    h1D_res_x_top_cut.append(  r.TH1F(laypfx + "h1D_res_x_top_cut", typePrefix+laypfx + "x residuals for top hit x;x residual (cm);counts",  100,-resRange,resRange)  )
    h1D_res_x_bottom_cut.append(  r.TH1F(laypfx + "h1D_res_x_bottom_cut", typePrefix+laypfx + "x residuals for bottom hit x;x residual (cm);counts",  100,-resRange,resRange)  )
#
    #h1D_res_x_actual_x.append(  r.TProfile(laypfx + "h1D_res_x_actual_x", typePrefix+laypfx + "x residual vs x;x (cm);x residual (cm)",  100,-200,200)  )
    h1D_res_x_actual_x.append(  r.TProfile(laypfx + "h1D_res_x_actual_x", typePrefix+laypfx + "x residual vs hit x;x (cm);x residual (cm)",  250,-XMax,XMax)  )
    h1D_res_x_track_x.append(  r.TProfile(laypfx + "h1D_res_x_track_x", typePrefix+laypfx + "x residual vs track x;x (cm);x residual (cm)",  250,-XMax,XMax)  )
    h1D_res_x_track_x_1st.append(  r.TProfile(laypfx + "h1D_res_x_track_x_1st", typePrefix+laypfx + "x residual vs track x;x (cm);x residual (cm)",  250,-XMax,XMax)  )
    h1D_res_x_track_x_2nd.append(  r.TProfile(laypfx + "h1D_res_x_track_x_2nd", typePrefix+laypfx + "x residual vs track x;x (cm);x residual (cm)",  250,-XMax,XMax)  )
    h1D_res_x_track_y.append(  r.TProfile(laypfx + "h1D_res_x_track_y", typePrefix+laypfx + "x residual vs track y;x (cm);x residual (cm)",  250,-XMax,XMax)  )
    h1D_pull_actual_x.append(  r.TProfile(laypfx + "h1D_pull_actual_x", typePrefix+laypfx + "pull of x residual vs hit x;x (cm);x residual (cm)",  100,-XMax,XMax)  )
    h1D_pull_tracks_x.append(  r.TProfile(laypfx + "h1D_pull_tracks_x", typePrefix+laypfx + "pull of x residual vs tracks x;x (cm);x residual (cm)",  100,-XMax,XMax)  )
    h1D_res_x_actual_x_TH2F.append(  r.TH2F(laypfx + "h1D_res_x_actual_x_TH2F", typePrefix+laypfx + "residual vs hit x;x (cm);x residual (cm)",  250,-XMax,XMax, 100, -10,10)  )
    h1D_res_x_track_x_TH2F.append(  r.TH2F(laypfx + "h1D_res_x_track_x_TH2F", typePrefix+laypfx + "residual vs track x;x (cm);x residual (cm)",  250,-XMax,XMax, 100, -10,10)  )
#    h1D_res_x_rproj.append(  r.TProfile(laypfx + "h1D_res_x_rproj", typePrefix+laypfx + "x residual vs r;r (cm);x residual (cm)",  26,-80.0,80.0)  )
#    #h1D_res_x_xproj.append(  r.TProfile(laypfx + "h1D_res_x_xproj", laypfx + "x residual vs x",  nXbins,-100.0,100.0)  )
#    # rphiproj has 50 bins, but rproj has 26 bins. This is so that the number of bins for PHYSICAL regions is the same
#    h1D_res_x_rphiproj.append(  r.TProfile(laypfx + "h1D_res_x_rphiproj", typePrefix+laypfx + "x residual vs scaled-phi;scaled-phi (cm);x residual (cm)",  50,-80.0,80.0)  )
#
#    # there are 2.530335 cm between each of the 48 wgs, so to get a binsize of 2.530335, we need (80-(-80))/2.530335=63.2=64 bins
#    h1D_res_rphi_yproj.append(  r.TProfile(laypfx + "h1D_res_rphi_yproj", typePrefix+laypfx + "rphi residual vs local y;local y (cm);rphi residual (cm)",  64,-80.0,80.0)  )
#
#    #h1D_res_x[i].StatOverflows(True)
#    #h1D_res_x_rproj[i].StatOverflows(True)
#    #h1D_res_x_rphiproj[i].StatOverflows(True)
#
    if(i < 4):
        h2D_res_y_actual.append(  r.TProfile2D(laypfx + "h2D_res_y_actual", typePrefix+laypfx + "avg y res at actual hit positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
        h2D_res_y_tracks.append(  r.TProfile2D(laypfx + "h2D_res_y_tracks", typePrefix+laypfx + "avg y res at tracks hit positions;local x (cm);local y (cm)",  nXbins,-XMax,XMax,  int(nXbins*xyRatio), -YMax,YMax)  )
        h1D_res_y.append(  r.TH1F(laypfx + "h1D_res_y", typePrefix+laypfx + "y residuals;y residual (cm);counts",  100,-resRange,resRange)  )

fh = r.TFile(filename)
tt = r.gDirectory.Get("dt_layer_ttree")
#
nPosMu = 0
nNegMu = 0
nMuons = 0
posIndexPt = []
negIndexPt = []
minPt = 9999.0


#for i,muon in enumerate(tt):
#    if( nMuons >= 100000): break
#    nMuons += 1
#    pt = muon.pz
#    if(not muon.select or muon.nlayers < 3): continue
#    wheel, station, sector = ord(muon.wheel), ord(muon.station), ord(muon.sector)
#    if (not (selectedDT == (wheel, station, sector))): continue
#
#            
#    if muon.charge == 1:
#        nPosMu += 1
#        posIndexPt.append([i,pt])
#    else:
#        nNegMu += 1
#        negIndexPt.append([i,pt])
#    if(pt < minPt): minPt = pt
#
#print nMuons
#goodIndices = getGoodMuonIndices(posIndexPt, negIndexPt, 17, minPt)
#print ">>> Found %i positive muons and %i negative muons" % (nPosMu, nNegMu)
#nMuonsToConsider = len(goodIndices) 
#print ">>> Considering ",nMuonsToConsider,"total muons after {q/pT,q/pz} equalization"
#
#goodMuons = {}
#for i in range(nMuons): goodMuons[i] = 0
#for e in goodIndices: goodMuons[e] = 1
TProfile_phi_view = r.TProfile("h1D_phi_view", "h1D_phi_view", 200, -3.2, 3.2, - 50, 50)
#
count = 0
for idx, muon in enumerate(tt):
    count += 1
    if(count % 2000 == 0): print ">>>",count 
#
    #if(count > 1000000): break
    if(count > count_max and count_max>0): break
    if(not muon.select or muon.nlayers < 6): continue

    wheel, station, sector = ord(muon.wheel), ord(muon.sector), ord(muon.station)

    res_x = muon.track_x[2]-muon.hit_x[2] 
    TProfile_phi_view.Fill(muon.phi, res_x)
    if(sector ==2):
        h2D_res_x_phi_eta.Fill(muon.phi, muon.eta, res_x)
        h2D_res_y_phi_eta.Fill(muon.phi, muon.eta, res_x)
    if(sector ==4 and muon.hit_x[0] >-500):
        hit = muon.hit_x[0]+station*450
        offset = 125
        if( station > 4):
            hit = hit+250

        if(station == 4 ): 
            hit = hit -offset + 75
            h1D_x_occupancy_whole_ring_special_chambers.Fill(hit)
            h2D_hit_x_res_x_whole_ring_special_chambers.Fill(hit, res_x)
        if(station == 13 ): 
            hit = muon.hit_x[0]+4*450  + 2.2*offset
            h1D_x_occupancy_whole_ring_special_chambers.Fill(hit)
            h2D_hit_x_res_x_whole_ring_special_chambers.Fill(hit, res_x)
        if(station ==10 ):
            hit = hit -offset
            h1D_x_occupancy_whole_ring_special_chambers.Fill(hit)
            h2D_hit_x_res_x_whole_ring_special_chambers.Fill(hit, res_x)
        if(station == 14 ): 
            hit = muon.hit_x[0]+10*475 + offset*1.1
            h1D_x_occupancy_whole_ring_special_chambers.Fill(hit)
            h2D_hit_x_res_x_whole_ring_special_chambers.Fill(hit, res_x)

        h2D_hit_x_res_x_whole_ring.Fill(hit,res_x)
        h1D_x_occupancy_whole_ring.Fill(hit)


    #if(selectedDT == (wheel, sector, station) and abs(muon.hit_x[4]) < 90):
    if(selectedDT == (wheel, sector, station) and muon.track_x[2] < fid_cut_x_1 and muon.track_x[2] > fid_cut_x_2 and muon.track_y_x_layer[2] < fid_cut_y_1 and muon.track_y_x_layer[2] > fid_cut_y_2):
    #if(selectedDT[1] == sector):
    
#        if(goodMuons[idx] != 1): continue
        #if(not muon.select): continue
        h1D_pt.Fill(muon.pt)
        h1D_pz.Fill(muon.pz)
        h1D_p.Fill(math.sqrt(muon.pt*muon.pt+muon.pz*muon.pz))
        h1D_eta.Fill(muon.eta)
        phi = muon.phi
        #if(muon.phi < 0): phi += 2*3.14159265
        h1D_phi.Fill(phi)
        for i in range(NUMLAYERS_X):
            try:
                actual_x, actual_y = muon.hit_x[i], muon.hit_y[i]
                res_x, res_y = muon.track_x[i]-muon.hit_x[i], muon.track_y_x_layer[i]-muon.hit_y[i]
                track_x = muon.track_x[i]
                track_y = muon.track_y_x_layer[i]
            except:
                layStr = "muon.lay%i_" % (i+1)
#                actual_x, actual_y = eval(layStr+"x"), eval(layStr+"y")
#                res_x, res_y = eval(layStr+"res_x"), eval(layStr+"res_y")
            #print i, res_x, track_x,  actual_x
            if(abs(track_x) >300): continue
            if(abs(actual_x) > 300): continue
            #if(abs(track_y) > 300): continue
            #if(abs(actual_y) > 300): continue
            # roughly discard tails
            if(abs(res_x) > resRange): continue 

            # propagated track = actual + residual

            #res_x = 11.0
            #track_x = 100.0
            #actual_x = track_x - res_x
            #adjust res_x for tilting bias
            
            adjusted_track = track_x -100
            adjusted_actual = actual_x -100
            radius_for_adjusted = math.sqrt(adjusted_track*adjusted_track+ DT_0_4_xpinR*DT_0_4_xpinR)
            #assume scatters at half distance:
            radius_for_adjusted_scatter = radius_for_adjusted/2

            phi_for_track = math.atan(adjusted_track/DT_0_4_xpinR)

            x1 = math.sin(phi_for_track)*radius_for_adjusted_scatter
            y1 = math.cos(phi_for_track)*radius_for_adjusted_scatter

            phi_for_hit = math.atan((actual_x-x1)/(DT_0_4_xpinR-y1))
            theta_for_adjusted = phi_for_hit - phi_for_track
            #adjusted_res_x = res_x*(1-phi_for_track*math.tan(theta_for_adjusted))
            #adjusted_res_x = -1.0/math.tan(theta_for_adjusted)*(radius_for_adjusted-radius_for_adjusted_scatter)
            #adjusted_res_x = res_x*(1-math.atan(track_x/DT_0_4_xpinR)*res_x/(math.sqrt(DT_0_4_xpinR*DT_0_4_xpinR + track_x*track_x)))
            adjusted_track_x = phi_for_track*DT_0_4_xpinR
            adjusted_hit_x = math.atan(adjusted_actual/DT_0_4_xpinR)*DT_0_4_xpinR   
            adjusted_res_x = adjusted_track_x - adjusted_hit_x
#            height_scatter = DT_0_4_xpinR-y1
#            hit_x_sctatter_frame = actual_x - x1
#            slope = (x1-actual_x)/(y1-DT_0_4_xpinR)
#            line_constant = DT_0_4_xpinR/(slope*actual_x)
#            adjusted_hit_x = 
            #print adjusted_res_x, "\t", res_x,"\t", theta_for_adjusted, "\t", radius_for_adjusted, "\t", track_x, "\t", actual_x
#            # do translations or rotations
#
#            #actual_x, actual_y = 0.99*actual_x, 0.99*actual_y # make chamber 1% smaller wrt to x,y=0,0
#            #actual_x, actual_y = actual_x, actual_y+0.30 # shift y up by 3 mm
#            #actual_x, actual_y = applyRotation(actual_x, actual_y, 0.001) # rotate CCW by 1 mrad
#            #actual_x, actual_y = actual_x + 0.03, actual_y # shift x up by 300 microns
#            # recalculate residuals
#            #res_x, res_y = track_x - actual_x, track_y - actual_y
#
#            ######
#
#            angle = math.atan(1.0 * actual_x / (ME13pinR + actual_y))
#            angle_track = math.atan(1.0 * track_x / (ME13pinR + track_y))
#            rphi_track = (ME13pinR)*math.atan(track_x / (ME13pinR + track_y))
#            res_rphi = math.cos(angle)*res_x + math.sin(angle)*res_y
#        
#            if(i==2): # if layer 3 ("center" of chamber)
#                h1D_actual_localy.Fill(actual_y)
#                h1D_tracks_localy.Fill(track_y)
#            
#                h1D_actual_angle.Fill(angle)
#                h1D_tracks_angle.Fill(angle_track)
#
#                h2D_nlayers_hit.Fill(actual_x, actual_y, muon.nlayers)
#                h2D_nDT_hit.Fill(actual_x, actual_y, muon.nDT)
#                h2D_nDT_hit.Fill(actual_x, actual_y, muon.nDT)
#                h2D_nDTDT_hit.Fill(actual_x, actual_y, muon.nDT+muon.nDT)
#                h2D_nTracker_hit.Fill(actual_x, actual_y, muon.nTracker)
#
#
            h1D_res_x[i].Fill(res_x)
            if(abs(track_x) < 90):
                h1D_res_x_center_cut[i].Fill(res_x)
            if((track_x) > 90):
                h1D_res_x_top_cut[i].Fill(res_x)
            if((track_x) < -90):
                h1D_res_x_bottom_cut[i].Fill(res_x)
            h1D_res_x_actual_x[i].Fill(actual_x,res_x) 
            h1D_res_x_track_x[i].Fill(track_x,res_x) 
            if(count%2 == 0):
                h1D_res_x_track_x_1st[i].Fill(track_x,res_x) 
            else:
                h1D_res_x_track_x_2nd[i].Fill(track_x,res_x) 

            h1D_res_x_track_y[i].Fill(track_y,res_x) 
            h1D_pull_tracks_x[i].Fill(track_x,res_x) 
            h1D_pull_actual_x[i].Fill(actual_x,res_x) 
            h1D_res_x_actual_x_TH2F[i].Fill(actual_x,res_x) 
            h1D_res_x_track_x_TH2F[i].Fill(track_x,res_x) 
#            h2D_cnt_actual[i].Fill(actual_x, actual_y)
#            h2D_cnt_tracks[i].Fill(track_x, track_y)
            h2D_res_x_actual[i].Fill(actual_x, actual_y, res_x)
            h2D_res_x_tracks[i].Fill(track_x, track_y, res_x)
            if(i < 4):
                h1D_res_y[i].Fill(res_y)
                h2D_res_y_actual[i].Fill(actual_x, actual_y, res_y)
                h2D_res_y_tracks[i].Fill(track_x, track_y, res_y)
            if(i==1):
                h2D_resx_res_y.Fill(res_x,res_y)
                h1D_x_occupancy.Fill(actual_x)
                #print int((actual_x + XMax)/(2*XMax/n_x_chamber_cut)), actual_x
                res_x_bin = int((track_x + XMax)/(2*XMax/n_x_chamber_cut))
                if(res_x_bin >0):
                    h1D_res_x_cuts[res_x_bin].Fill(res_x)
                    h1D_adjusted_res_x_cuts[res_x_bin].Fill(adjusted_res_x)
                h1D_res_x_cuts[n_x_chamber_cut].Fill(res_x)
                h1D_adjusted_res_x_cuts[n_x_chamber_cut].Fill(adjusted_res_x)
                h1D_res_x_cuts_sym[n_x_chamber_cut].Fill(-res_x)
                h1D_adjusted_res_x_cuts_sym[n_x_chamber_cut].Fill(-adjusted_res_x)

                #print adjusted_res_x - res_x

#            h2D_pull_tracks[i].Fill(track_x, track_y, res_x)
#
#            h1D_res_x_rproj[i].Fill(rProjection(actual_y,angle), res_x)
#            h1D_res_x_rphiproj[i].Fill(rphi_track, res_x)
#            h1D_res_rphi_yproj[i].Fill(actual_y, res_rphi)
       

c1.SetRightMargin(0.32);
c1.SetGridx()
c1.SetGridy()
c1.SetCanvasSize(696,472)


os.system("mkdir -p " + prefix)
os.system("mkdir -p " + prefix + "/res_x_distributions_v_hit_x")
os.system("cp -p " + " indexbase.php " + prefix+"_index.php") #
os.system("cp -p " + " indexbase.php " + prefix+"/res_x_distributions_v_hit_x/"+"_index.php") #

r.gStyle.SetOptStat("rme")

#layerRotationR = []
#layerTranslation = []

#for i in range(NUMLAYERS_X):
#for i in range(NUMLAYERS_X):
for i in range(1):

    r.gStyle.SetPalette(1)
    r.gStyle.SetOptFit(0)
    #suffix = "_LAY" + str(i+1) + ".C"
    suffix = "_LAY" + str(i+1) + ".png"
    #suffix =  ".png"
    
    if(i == 0):
        h1D_pt.Draw()
        c1.SaveAs(prefix + "h1D_pt" + suffix)

        
        h1D_pz.Draw()
        c1.SaveAs(prefix + "h1D_pz" + suffix)

        h1D_p.Draw()
        c1.SaveAs(prefix + "h1D_p" + suffix)

        h1D_phi.Draw()
        zoomXrange(h1D_phi)
        c1.SaveAs(prefix + "h1D_phi" + suffix)

        zoomXrange(h1D_eta)
        h1D_eta.Draw()
        c1.SaveAs(prefix + "h1D_eta" + suffix)


        c1.SetCanvasSize(1000,672)

        TProfile_phi_view.Draw()
        c1.SaveAs(prefix + "TProfile_phi_view" + suffix)

        c1.SetCanvasSize(555,672)
        
        c1.SetCanvasSize(696,472)


    c1.SetCanvasSize(555,672)

    
    levels = [0.00, 0.50, 1.00]
    red    = [0.00, 0.50, 1.00]
    green  = [0.00, 0.50, 0.00]
    blue   = [1.00, 0.50, 0.00]
    #levels = [0.00, 0.35, 0.58, 0.78, 1.00]
    #red    = [0.00, 0.00, 0.87, 1.00, 0.51]
    #green  = [0.00, 0.00, 0.00, 0.00, 0.00]
    #blue   = [0.51, 1.00, 0.12, 0.00, 0.00]
    levels = array('d', levels)
    ncontours = 100

    r.TColor.CreateGradientColorTable(len(levels), levels, array('d', red), array('d', green), array('d', blue), ncontours)
    r.gStyle.SetNumberContours(ncontours)

    if(i == 0):
    

        c1.SetCanvasSize(2000,1000)
        h2D_res_y_phi_eta.GetZaxis().SetRangeUser(-3.5, 3.5)
        fixFlooredBins(h2D_res_y_phi_eta, minZ=-3.5)
        h2D_res_y_phi_eta.Draw("colz")
        c1.SaveAs(prefix + "h2D_res_y_phi_eta" + suffix)

        h2D_res_x_phi_eta.GetZaxis().SetRangeUser(-3.5, 3.5)
        fixFlooredBins(h2D_res_x_phi_eta, minZ=-3.5)
        h2D_res_x_phi_eta.Draw("colz")
        c1.SaveAs(prefix + "h2D_res_x_phi_eta" + suffix)
    
        leg = r.TLegend(0.1,0.85,0.2,0.9)
        h1D_x_occupancy_whole_ring_special_chambers.SetLineColor(2)
        h1D_x_occupancy_whole_ring.Draw()
        h1D_x_occupancy_whole_ring_special_chambers.Draw("same")
        leg.AddEntry(h1D_x_occupancy_whole_ring_special_chambers, "4,13,10 and 14")
        leg.Draw()
        c1.SetCanvasSize(4000,1000)
        c1.SaveAs(prefix + "h1D_x_occupancy_whole_ring" + suffix)


        h2D_hit_x_res_x_whole_ring_special_chambers.SetLineColor(2)
        h2D_hit_x_res_x_whole_ring.Draw("colz")
        #h2D_hit_x_res_x_whole_ring_special_chambers.Draw("same")
        leg.Draw()
        c1.SaveAs(prefix + "h2D_hit_x_res_x_whole_ring" + suffix)



        c1.SetCanvasSize(696,472)
        h2D_resx_res_y.Draw("colz")
        c1.SaveAs(prefix + "h2D_resx_res_y" + suffix)

        h1D_x_occupancy.Draw()
        c1.SaveAs(prefix + "h1D_x_occupancy" + suffix)



    h2D_res_x_actual[i].GetZaxis().SetRangeUser(-3.5, 3.5)
    fixFlooredBins(h2D_res_x_actual[i], minZ=-3.5)
    h2D_res_x_actual[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_res_x_actual" + suffix)

    h2D_res_x_tracks[i].GetZaxis().SetRangeUser(-1, 1)
    fixFlooredBins(h2D_res_x_tracks[i], minZ=-3.5)
    h2D_res_x_tracks[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_res_x_tracks" + suffix)
    h2D_res_x_tracks[i].Draw("CONTZ")
    c1.SaveAs(prefix + "h2D_res_x_tracks_contz" + suffix)
    h2D_res_x_tracks[i].Rebin2D(7,7)
    h2D_res_x_tracks[i].Draw("colz")
    c1.SaveAs(prefix + "h2D_res_x_tracks_small" + suffix)


    c1.SetCanvasSize(696,472)


    r.gStyle.SetOptFit(1) # display fitting parameters
# 
    h1D_res_x[i].Draw()
    fitCut(h1D_res_x[i], nSigma, "QC")
    c1.SaveAs(prefix + "h1D_res_x" + suffix)

    

    #for x in range(h2D_res_x_mean_v_hit_x.GetNbinsX()):
    #    h2D_res_x_mean_v_hit_x.Fill(x,x)
    #for x in range(h2D_res_x_mean_v_hit_x.GetNbinsX()):
    #    h2D_res_x_mean_v_hit_x.SetBinContent(x,0)

    #c1.SetLogy(1)
    if(i == 0):
        for j in range(n_x_chamber_cut):
            h1D_res_x_cuts[j].Draw()
            title = "/res_x_distributions_v_hit_x/h1D_res_x_cuts_%d_%d" % (j*2*XMax/n_x_chamber_cut-XMax, (j+1)*2*XMax/n_x_chamber_cut-XMax)
            #print title
            fitCut(h1D_res_x_cuts[j], nSigma, "QC")
            c1.SaveAs(prefix + title  + suffix)
            #print h1D_res_x_cuts[j].GetMean(), h1D_res_x_cuts[j].GetMeanError(), (j*2*XMax/n_x_chamber_cut-XMax+ (j+1)*2*XMax/n_x_chamber_cut-XMax)/2
            x_pos = (j*2*XMax/n_x_chamber_cut-XMax+ (j+1)*2*XMax/n_x_chamber_cut-XMax)/2
            bin_position = int(tProfile_cut_bins*(x_pos+XMax)/(2*XMax))
            #print "mean: " , h1D_res_x_cuts[j].GetMean()
           # print "skewness: " , h1D_res_x_cuts[j].GetSkewness()
           # print "RMS: " , h1D_res_x_cuts[j].GetRMS()


            if(j < n_x_chamber_cut):
                #c1.SetLogy(0)
                
                temp_histo =  h1D_adjusted_res_x_cuts[j].Clone()
                temp_histo.Divide(h1D_adjusted_res_x_cuts[j+1])
                temp_histo.GetXaxis().SetRangeUser(-4,4)
                temp_histo.Draw("E")
                #fitCut(h1D_adjusted_res_x_cuts[j], nSigma, "QC")
                title = "/res_x_distributions_v_hit_x/h1D_adjusted_res_x_cuts_divied_%d_%d" % (j*2*XMax/n_x_chamber_cut-XMax, (j+1)*2*XMax/n_x_chamber_cut-XMax)
                c1.SaveAs(prefix + title  + suffix)
                temp_histo.Reset()
                #c1.SetLogy(1)
            h1D_adjusted_res_x_cuts[j].Draw("E")
            fitCut(h1D_adjusted_res_x_cuts[j], nSigma, "QC")
            title = "/res_x_distributions_v_hit_x/h1D_adjusted_res_x_cuts_%d_%d" % (j*2*XMax/n_x_chamber_cut-XMax, (j+1)*2*XMax/n_x_chamber_cut-XMax)
            c1.SaveAs(prefix + title  + suffix)


            #fitparamsGauss = getFitParamsGauss(h1D_res_x_cuts[j])
            

            if(h1D_res_x_cuts[j].GetEntries() >300):

                h2D_adjusted_res_x_mean_v_hit_x.SetBinContent(bin_position, h1D_adjusted_res_x_cuts[j].GetMean())
                h2D_adjusted_res_x_mean_v_hit_x.SetBinError(bin_position, h1D_adjusted_res_x_cuts[j].GetMeanError())

                nq = 4
                xq = array('d', [0.] * nq)
                yq1 = array('d', [0.] * nq)
                yq2 = array('d', [0.] * nq)
                for i in xrange(nq):
                    xq[i] = float(i + 1) / nq
                
                h1D_adjusted_res_x_cuts[j].GetQuantiles(nq, yq2, xq)
                h2D_adjusted_res_x_median_v_hit_x.SetBinContent(bin_position, yq2[1])
                #error for median:http://davidmlane.com/hyperstat/A106993.html
                h2D_adjusted_res_x_median_v_hit_x.SetBinError(bin_position, h1D_adjusted_res_x_cuts[j].GetMeanError()*1.253)

                fitparamsGauss = getFitParamsGauss(h1D_adjusted_res_x_cuts[j])
                h2D_adjusted_res_x_fit_p0_v_hit_x.SetBinContent(bin_position, fitparamsGauss[0])
                h2D_adjusted_res_x_fit_p0_v_hit_x.SetBinError(bin_position, fitparamsGauss[1])

                h2D_adjusted_res_x_fit_p1_v_hit_x.SetBinContent(bin_position, fitparamsGauss[2])
                h2D_adjusted_res_x_fit_p1_v_hit_x.SetBinError(bin_position, fitparamsGauss[3])

                skewness = (h1D_adjusted_res_x_cuts[j].GetSkewness(1), h1D_adjusted_res_x_cuts[j].GetSkewness(11))
                h2D_adjusted_res_x_skewness_v_hit_x.SetBinContent(bin_position, skewness[0])
                h2D_adjusted_res_x_skewness_v_hit_x.SetBinError(bin_position, skewness[1])

                h2D_adjusted_res_x_RMS_v_hit_x.SetBinContent(bin_position, h1D_adjusted_res_x_cuts[j].GetRMS())
                h2D_adjusted_res_x_RMS_v_hit_x.SetBinError(bin_position, h1D_adjusted_res_x_cuts[j].GetRMSError())

                h2D_res_x_mean_v_hit_x.SetBinContent(bin_position, h1D_res_x_cuts[j].GetMean())
                h2D_res_x_mean_v_hit_x.SetBinError(bin_position, h1D_res_x_cuts[j].GetMeanError())

                
                h1D_res_x_cuts[j].GetQuantiles(nq, yq1, xq)
                h2D_res_x_median_v_hit_x.SetBinContent(bin_position, yq1[1])
                h2D_res_x_median_v_hit_x.SetBinError(bin_position, h1D_res_x_cuts[j].GetMeanError()*1.253)

                fitparamsGauss = getFitParamsGauss(h1D_res_x_cuts[j])
                h2D_res_x_fit_p0_v_hit_x.SetBinContent(bin_position, fitparamsGauss[0])
                h2D_res_x_fit_p0_v_hit_x.SetBinError(bin_position, fitparamsGauss[1])

                h2D_res_x_fit_p1_v_hit_x.SetBinContent(bin_position, fitparamsGauss[2])
                h2D_res_x_fit_p1_v_hit_x.SetBinError(bin_position, fitparamsGauss[3])

                skewness = (h1D_res_x_cuts[j].GetSkewness(1), h1D_res_x_cuts[j].GetSkewness(11))
                h2D_res_x_skewness_v_hit_x.SetBinContent(bin_position, skewness[0])
                h2D_res_x_skewness_v_hit_x.SetBinError(bin_position, skewness[1])

                h2D_res_x_RMS_v_hit_x.SetBinContent(bin_position, h1D_res_x_cuts[j].GetRMS())
                h2D_res_x_RMS_v_hit_x.SetBinError(bin_position, h1D_res_x_cuts[j].GetRMSError())

            #h2D_res_x_skewness_v_hit_x.SetBinContent(bin_position, h1D_res_x_cuts[j].GetSkewness(1))
            #h2D_res_x_skewness_v_hit_x.SetBinError(bin_position, h1D_res_x_cuts[j].GetSkewness(11))
        #c1.SetLogy(0)

        print "skew res_x\t",  h1D_res_x_cuts[j+1].GetSkewness(1), "\t", h1D_res_x_cuts[j+1].GetSkewness(11), "\tskew adjusted res_x\t", h1D_adjusted_res_x_cuts[j+1].GetSkewness(1), "\t", h1D_adjusted_res_x_cuts[j+1].GetSkewness(11)

        temp_histo = h1D_res_x_cuts[j+1].Clone() 
        temp_histo.Divide(h1D_res_x_cuts_sym[j+1])
        temp_histo.Fit("pol3")
        temp_histo.Draw()
        title = "/res_x_distributions_v_hit_x/h1D_res_x_cuts_total_assymetry"
        c1.SaveAs(prefix + title  + suffix)
        temp_histo.Rebin(3)
        temp_histo.Scale(1.0/3)
        temp_histo.GetXaxis().SetRangeUser(-4,4)
        temp_histo.Fit("pol1")
        temp_histo.Draw()
        title = "/res_x_distributions_v_hit_x/h1D_res_x_cuts_total_assymetry_reduced_range"
        c1.SaveAs(prefix + title  + suffix)
            
        h1D_res_x_cuts[j+1].Draw()
        fitCut(h1D_res_x_cuts[j+1], nSigma, "QC")
        title = "/res_x_distributions_v_hit_x/h1D_res_x_cuts_total"
        c1.SaveAs(prefix + title  + suffix)

        temp_histo = h1D_adjusted_res_x_cuts[j+1].Clone() 
        temp_histo.Divide(h1D_adjusted_res_x_cuts_sym[j+1])
        temp_histo.Fit("pol3")
        temp_histo.Draw()
        title = "/res_x_distributions_v_hit_x/h1D_adjusted_res_x_cuts_total_assymetry"
        c1.SaveAs(prefix + title  + suffix)
        temp_histo.Rebin(3)
        temp_histo.Scale(1.0/3)
        temp_histo.GetXaxis().SetRangeUser(-4,4)
        temp_histo.Fit("pol1")
        temp_histo.Draw()
        title = "/res_x_distributions_v_hit_x/h1D_adjusted_res_x_cuts_total_assymetry_reduced_range"
        c1.SaveAs(prefix + title  + suffix)

        h1D_adjusted_res_x_cuts[j+1].Draw()
        fitCut(h1D_adjusted_res_x_cuts[j+1], nSigma, "QC")
        title = "/res_x_distributions_v_hit_x/h1D_adjusted_res_x_cuts_total"
        c1.SaveAs(prefix + title  + suffix)

        h1D_res_x_cuts[j+1].Draw()
        fitCut(h1D_res_x_cuts[j+1], nSigma, "QC")
        title = "/res_x_distributions_v_hit_x/h1D_res_x_cuts_total"
        c1.SaveAs(prefix + title  + suffix)

        print "res x skew\t", 
        h1D_adjusted_res_x_cuts[j+1].Draw()
        fitCut(h1D_adjusted_res_x_cuts[j+1], nSigma, "QC")
        title = "/res_x_distributions_v_hit_x/h1D_adjusted_res_x_cuts_total"
        c1.SaveAs(prefix + title  + suffix)

        title = "/res_x_distributions_v_hit_x/" 

        leg_adjusted = r.TLegend(0.1,0.85,0.2,0.9)

        h2D_adjusted_res_x_fit_p0_v_hit_x.SetLineColor(2)
        h2D_adjusted_res_x_fit_p1_v_hit_x.SetLineColor(2)
        h2D_adjusted_res_x_RMS_v_hit_x.SetLineColor(2)
        h2D_adjusted_res_x_skewness_v_hit_x.SetLineColor(2)
        h2D_adjusted_res_x_mean_v_hit_x.SetLineColor(2)
        h2D_adjusted_res_x_median_v_hit_x.SetLineColor(2)

        leg_adjusted.AddEntry(h2D_adjusted_res_x_fit_p0_v_hit_x, "adjusted data")
        leg_adjusted.AddEntry(h2D_res_x_fit_p0_v_hit_x, "original data")

        h2D_res_x_fit_p0_v_hit_x.Draw()
        c1.SaveAs(prefix + title + "h2D_res_x_fit_p0_v_hit_x" + suffix)
        h2D_adjusted_res_x_fit_p0_v_hit_x.Draw("same")
        leg_adjusted.Draw()
        c1.SaveAs(prefix + title + "h2D_adjusted_res_x_fit_p0_v_hit_x" + suffix)

        h2D_res_x_fit_p1_v_hit_x.SetMaximum(.08)
        h2D_res_x_fit_p1_v_hit_x.SetMinimum(-.08)
        h2D_res_x_fit_p1_v_hit_x.Fit("pol1")
        h2D_res_x_fit_p1_v_hit_x.Draw()
        c1.SaveAs(prefix + title + "h2D_res_x_fit_p1_v_hit_x" + suffix)
        h2D_adjusted_res_x_fit_p1_v_hit_x.Draw("same")
        leg_adjusted.Draw()
        c1.SaveAs(prefix + title + "h2D_adjusted_res_x_fit_p1_v_hit_x" + suffix)

        h2D_res_x_RMS_v_hit_x.Draw()
        c1.SaveAs(prefix + title + "h2D_res_x_RMS_v_hit_x" + suffix)
        h2D_adjusted_res_x_RMS_v_hit_x.Draw("same")
        leg_adjusted.Draw()
        c1.SaveAs(prefix + title + "h2D_adjusted_res_x_RMS_v_hit_x" + suffix)

        h2D_res_x_skewness_v_hit_x.Fit("pol1")
        h2D_res_x_skewness_v_hit_x.Draw()
        c1.SaveAs(prefix + title + "h2D_res_x_skewness_v_hit_x" + suffix)
        h2D_adjusted_res_x_skewness_v_hit_x.Draw("same")
        leg_adjusted.Draw()
        c1.SaveAs(prefix + title + "h2D_adjusted_res_x_skewness_v_hit_x" + suffix)

        h2D_res_x_median_v_hit_x.SetMaximum(.08)
        h2D_res_x_median_v_hit_x.SetMinimum(-.08)
        h2D_res_x_median_v_hit_x.Fit("pol1")
        h2D_res_x_median_v_hit_x.Draw()
        c1.SaveAs(prefix + title + "h2D_res_x_median_v_hit_x" + suffix)
        h2D_adjusted_res_x_median_v_hit_x.Draw("same")
        leg_adjusted.Draw()
        c1.SaveAs(prefix + title + "h2D_adjusted_res_x_median_v_hit_x" + suffix)

        h2D_res_x_mean_v_hit_x.SetMaximum(.08)
        h2D_res_x_mean_v_hit_x.SetMinimum(-.08)
        h2D_res_x_mean_v_hit_x.Fit("pol1")
        h2D_res_x_mean_v_hit_x.Draw()
        c1.SaveAs(prefix + title + "h2D_res_x_mean_v_hit_x" + suffix)
        h2D_adjusted_res_x_mean_v_hit_x.Draw("same")
        leg_adjusted.Draw()
        c1.SaveAs(prefix + title + "h2D_adjusted_res_x_mean_v_hit_x" + suffix)
    
    h1D_res_x_center_cut[i].Draw()
    fitCut(h1D_res_x_center_cut[i], nSigma, "QC")
    c1.SaveAs(prefix + "h1D_res_x_center_cut" + suffix)

    h1D_res_x_top_cut[i].Draw()
    fitCut(h1D_res_x_top_cut[i], nSigma, "QC")
    c1.SaveAs(prefix + "h1D_res_x_top_cut" + suffix)

    h1D_res_x_bottom_cut[i].Draw()
    fitCut(h1D_res_x_bottom_cut[i], nSigma, "QC")
    c1.SaveAs(prefix + "h1D_res_x_bottom_cut" + suffix)

    #h1D_pull_tracks_x[i].Rebin(5)
    #h1D_pull_actual_x[i].Rebin(5)
    for x in range(h1D_pull_tracks_x[i].GetNbinsX()):
        mu = h1D_pull_tracks_x[i].GetBinContent(x)
        entries = h1D_pull_tracks_x[i].GetBinEntries(x)
        if(entries > 0): h1D_pull_tracks_x[i].SetBinContent(x,mu*entries*entries)
 #       print "pull at bin x:\t" , x, "\t", mu, "\t", entries, "\t", mu*entries*entries
        mu = h1D_pull_actual_x[i].GetBinContent(x)
        entries = h1D_pull_actual_x[i].GetBinEntries(x)
        if(entries > 0): h1D_pull_actual_x[i].SetBinContent(x,mu*entries*entries)



    c1.SetCanvasSize(2000,1000)
    h1D_pull_tracks_x[i].Fit("pol2")
    h1D_pull_tracks_x[i].Draw()
    c1.SaveAs(prefix + "h1D_pull_tracks_x" + suffix)

    h1D_pull_actual_x[i].Draw()
    h1D_pull_actual_x[i].Fit("pol2")
    c1.SaveAs(prefix + "h1D_pull_actual_x" + suffix)


    h1D_res_x_actual_x[i].Rebin(10)
    h1D_res_x_track_x[i].Rebin(10)
    h1D_res_x_track_x[i].SetMinimum(-.08)
    h1D_res_x_track_x[i].SetMaximum(.08)
    h1D_res_x_actual_x[i].Draw()
    c1.SaveAs(prefix + "h1D_res_x_actual_x" + suffix)
    h1D_res_x_track_x[i].Draw()
    c1.SaveAs(prefix + "h1D_res_x_track_x" + suffix)

    h1D_res_x_track_x_1st[i].Rebin(10)
    h1D_res_x_track_x_1st[i].SetLineColor(2)
    h1D_res_x_track_x_1st[i].SetMinimum(-.08)
    h1D_res_x_track_x_1st[i].SetMaximum(.08)
    h1D_res_x_track_x_1st[i].Draw()
    c1.SaveAs(prefix + "h1D_res_x_track_x_1st" + suffix)

    h1D_res_x_track_x_2nd[i].Rebin(10)
    h1D_res_x_track_x_2nd[i].SetMinimum(-.08)
    h1D_res_x_track_x_2nd[i].SetMaximum(.08)
    h1D_res_x_track_x_2nd[i].Draw()
    c1.SaveAs(prefix + "h1D_res_x_track_x_2nd" + suffix)
    h1D_res_x_track_x_2nd[i].Draw()
    h1D_res_x_track_x_1st[i].Draw("same")
    c1.SaveAs(prefix + "h1D_res_x_track_x_both" + suffix)

    h1D_res_x_track_y[i].Rebin(10)
    h1D_res_x_track_y[i].Draw()
    c1.SaveAs(prefix + "h1D_res_x_track_y" + suffix)

    h1D_res_x_actual_x_TH2F[i].Draw("colz")
    c1.SaveAs(prefix + "h1D_res_x_actual_x_TH2F" + suffix)
    h1D_res_x_actual_x_TH2F[i].Draw("CANDLEX1")
    h1D_res_x_actual_x_TH2F[i].Rebin2D(3,1) 
    c1.SaveAs(prefix + "h1D_res_x_actual_x_TH2F_candle" + suffix)
    h1D_res_x_actual_x_TH2F[i].Rebin2D(1,3) 
    h1D_res_x_actual_x_TH2F[i].Draw("CONT1")
    c1.SaveAs(prefix + "h1D_res_x_actual_x_TH2F_cont" + suffix)
    h1D_res_x_actual_x_TH2F[i].Rebin2D(4,1) 
    h1D_res_x_actual_x_TH2F[i].Draw("CANDLEX1")
    c1.SaveAs(prefix + "h1D_res_x_actual_x_TH2F_candle2" + suffix)
    h1D_res_x_track_x_TH2F[i].Draw("colz")
    c1.SaveAs(prefix + "h1D_res_x_track_x_TH2F_colz" + suffix)
    h1D_res_x_track_x_TH2F[i].Rebin2D(15,1)
    h1D_res_x_track_x_TH2F[i].Draw("CANDLEX1")
    c1.SaveAs(prefix + "h1D_res_x_track_x_TH2F_candle2" + suffix)
    c1.SetCanvasSize(696,472)

    if(i<4):

        h2D_res_y_actual[i].GetZaxis().SetRangeUser(-3.5, 3.5)
        fixFlooredBins(h2D_res_y_actual[i], minZ=-3.5)
        h2D_res_y_actual[i].Draw("colz")
        c1.SaveAs(prefix + "h2D_res_y_actual" + suffix)

        h2D_res_y_tracks[i].GetZaxis().SetRangeUser(-3.5, 3.5)
        fixFlooredBins(h2D_res_y_tracks[i], minZ=-3.5)
        h2D_res_y_tracks[i].Draw("colz")
        c1.SaveAs(prefix + "h2D_res_y_tracks" + suffix)

        h1D_res_y[i].Draw()
        fitCut(h1D_res_y[i], nSigma, "QC")
        c1.SaveAs(prefix + "h1D_res_y" + suffix)

r.gStyle.SetPalette(1)

