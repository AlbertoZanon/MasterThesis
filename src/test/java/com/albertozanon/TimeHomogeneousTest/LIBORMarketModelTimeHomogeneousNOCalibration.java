/*
 * (c) Copyright Christian P. Fries, Germany. Contact: email@christian-fries.de.
 *
 * Created on 16.01.2015
 */
package com.albertozanon.TimeHomogeneousTest;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.time.LocalDate;
import java.time.Month;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.FutureTask;

import org.junit.Assert;

import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.marketdata.calibration.ParameterObject;
import net.finmath.marketdata.calibration.Solver;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.AnalyticModelFromCurvesAndVols;
import net.finmath.marketdata.model.curves.Curve;
import net.finmath.marketdata.model.curves.CurveInterpolation.ExtrapolationMethod;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationEntity;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationMethod;
import net.finmath.marketdata.model.volatilities.AbstractVolatilitySurface;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterpolation;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveFromDiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.marketdata.products.AnalyticProduct;
import net.finmath.marketdata.products.Swap;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromTermStructureModel;
import net.finmath.montecarlo.interestrate.TermStructureModel;
import net.finmath.montecarlo.interestrate.models.LIBORMarketModelFromCovarianceModel;

import com.albertozanon.TimeHomogeneouLMM.CapletOnCompoundedBackward;
import com.albertozanon.TimeHomogeneouLMM.LIBORMarketModelWithTenorRefinement;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.BlendedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.DisplacedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelExponentialForm5Param;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelFourParameterExponentialForm;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelPiecewiseConstant;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelTimeHomogenousPiecewiseConstant;
import net.finmath.montecarlo.interestrate.models.covariance.TermStructCovarianceModelFromLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.TermStructureCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.TermStructureTenorTimeScalingInterface;
//import net.finmath.montecarlo.interestrate.models.covariance.TermStructureTenorTimeScalingPicewiseConstant;

import net.finmath.montecarlo.interestrate.models.covariance.TermStructureTenorTimeScalingPicewiseConstant;

import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.Caplet;
import net.finmath.montecarlo.interestrate.products.SwaptionSimple;
import net.finmath.montecarlo.interestrate.products.Caplet.ValueUnit;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.Optimizer;
import net.finmath.optimizer.OptimizerFactory;
import net.finmath.optimizer.OptimizerFactoryLevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.Schedule;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.time.TimeDiscretizationFromArray.ShortPeriodLocation;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;
import net.finmath.time.daycount.DayCountConvention_ACT_365;

/**
 * This class tests the LIBOR market model and products.
 *
 * @author Christian Fries
 */
public class LIBORMarketModelTimeHomogeneousNOCalibration {

	private final int numberOfPaths		= 2000;
	private final int numberOfFactors	= 1;

	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.0000%;-##0.0000%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterParam		= new DecimalFormat(" #0.00000; -#0.00000", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));

	
	public static void main(final String[] args) throws CalculationException, SolverException {
		final LIBORMarketModelTimeHomogeneousNOCalibration test = new LIBORMarketModelTimeHomogeneousNOCalibration();
		test.testATMSwaptionCalibration();
	}

	//constructor
	public LIBORMarketModelTimeHomogeneousNOCalibration() {}

		public void testATMSwaptionCalibration() throws CalculationException, SolverException {
		/*
		 * Calibration test
		 */
		System.out.println("Calibration to Swaptions:");

		final AnalyticModel curveModel = getCalibratedCurve();

		// Create the forward curve (initial value of the LIBOR market model)
		final ForwardCurve forwardCurve = curveModel.getForwardCurve("ForwardCurveFromDiscountCurve(discountCurve-EUR,1D)");

		final DiscountCurve discountCurve = curveModel.getDiscountCurve("discountCurve-EUR");

			double[] paramTimeScalingVolParam = { 0.09, 0.09, 0.09, 0.09, -0.001900851034611, -0.00419575274955029, -0.00442084917971079, -0.00424776052425492, 0.0013764351846531, -0.00424098708715426, -0.00432276020850956, -0.00484327828049977, -0.00350582260076357, 0.004103522850748, -0.0033822564678249, -0.00223657909794653, 0.00143883962989193, 0.00156242850096184, 0.000881836104654567, 0.000878353505005478, 0.00202916721414624, 0.00744677974748718, -0.00499777720394327, -0.0032375775025563, -0.00609439733504146, -0.00238708478822353, -0.00183051148761336, 0.00745903406524448, 0.00145028601939714, 0.00149330089465308, 0.00156737975294966, 0.00161404363808074, 0.00141793953028837, -0.0052268458202927, 0.00365425462942354, 0.00946325638395848, 0.00681170414948539, -0.000940158807611766, -0.00699719094804223, 0.0135496307619622, 0.0185688809967169, 0.0311470500521249, 0.0311443375311536, 0.0295441217680911, 0.011892840783714, 0.0105499541420107, 0.00980427909840784, 0.0108669835180354, -0.00103425214317539, -0.00132711749747869, -0.00136280556843658, -0.000492166651541481, -0.00684797314217561, -0.00716493781129585, -0.00716738015141388, -0.00795060850560716, -0.00805094521185637, -0.00798743825874311, -0.00774609653300004, -0.00857869815592934, 0.00482854937626569, 0.00197785983094334, 0.00280251478929372, 0.000213109666955375, -0.000596924822051647, 0.0123793060132364, 0.0157545703846876, -0.0011865071259416, -0.00559066498178182, -0.00458095318497953, -0.00532061024929078, -0.00671344086117528, -0.00900000000000005, -0.00721630270143222, -0.0058148332445711, -0.0065798467479648, -0.00833700942159353, -0.00478590931714208, -0.00679184440698776, -0.00865227523856503, -0.001362554514414, -0.00215421678503531, 0.00563395174074059, 0.00288830137075635, -0.00322590555903346, -0.00319071778582056, -0.00332317957646992, -0.00333655953685308, -0.00737945746580621, -0.00559783201459666, -0.0064796064969056, -0.00681079187296774, -0.00562770110641906, -0.00386385774918665, -0.00543514188911274, -0.00313767056339287, -0.00610601314773561, -0.00804130953605721, -0.0059856471488905, -0.0051825868236844, -0.00157742038394218, 0.00234885673931785, -0.00900000000000005, -0.00900000000000005, -0.00733946049440504, -0.00769440762248876, -0.00900000000000005, -0.00899617654774175, -0.00701160001265094, -0.00858582578346101, -0.00746246514379436, -0.00900000000000005, -0.00900000000000005, 0.000598322579645014, -0.00900000000000005, -0.00900000000000005, -0.00892773221649463, -0.00893756172546233, -0.00893395964777937, -0.00900000000000005, -0.00900000000000005, 0.00329899008574813, 0.00281093577300424, 0.00226508346888493, -0.00509846886780593, -0.00509846886780252, -0.00508223132113386, -0.00508223132113414, -0.00548312892577001, -0.0054831289257163, -0.0054950348199361, -0.00549503481997561, 0.0034570427157422, 0.00345704271559782, 0.00361842409308451, 0.00361842409316437, -0.0028284582869719, -0.0028284582867792, -0.00242531833401301, -0.00242531833380979, -0.00813797926683094, -0.00813797926681104, -0.00813699944417436, -0.00813699944416896, -0.00883451143820678, -0.00883451143820678, -0.0088348640358015, -0.00883486403580235, -0.00900000000000005, -0.00900000000000005, -0.00900000000000005, -0.00900000000000005, -0.00864709861410944, -0.00864709861410944, -0.00864798790067908, -0.00864798790067908, -0.00900000000000005, -0.00900000000000005, -0.00900000000000005, -0.00900000000000005, 0.004, -0.494432859305378, -0.494432859006778, -0.494432859307916, -0.00157456511415939, 1.21851090326318E-06, 0.0027202838940817, 0.00452697788832656, 0.00549810253240709, 0.010995027658776, 0.00333689489212563, 0.00197221982391949, 0.00141531726062687, -0.00188531949594071, 0.004, 0.004, 0.004, 0.004, -0.000141069120803055, 0.004, 0.004, 0.004, 0.004, -0.000387912613996856, 0.004, 0.004, 0.004, 0.004, -0.000270490730546011, 0.004, 0.004, 0.004, 0.004, -0.000938307181808428, 0.004, 0.004, 0.004, 0.004, -0.00113333707074173, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 1.10645690745991, 1.10645690723855, -0.000196522038430659, -0.000155975866775778, 0.00392811095533235, 0.00734455193779385, 0.00860604548768743, 0.0115863782454146, 0.0097501847215871, 0.0127380105057194, 0.0153888755271697, 0.000730431658179325, 0.004, 0.004, 0.004, 0.004, -0.000209317243014865, 0.004, 0.004, 0.004, 0.004, -0.000225865837606041, 0.004, 0.004, 0.004, 0.004, -0.000314735418794094, 0.004, 0.004, 0.004, 0.004, -0.000546643602310953, 0.004, 0.004, 0.004, 0.004, -0.00437598358073981, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.00389706087474156, 0.00375800084091632, 0.00384260482602033, 0.00487086422437719, 0.00440649505996209, 0.00326807613630755, 0.00363821688177974, 0.00178214622005669, 0.00459380269628534, 0.00337814431416731, 0.00148870087545478, 0.00347388051058275, -0.000150015736133534, 0.004, 0.004, 0.004, 0.00174036886717814, 0.004, 0.004, 0.004, 0.004, 0.00184685842633084, 0.004, 0.004, 0.004, 0.004, 0.00489521607373615, 0.004, 0.004, 0.004, 0.004, 0.000200747953398207, 0.004, 0.004, 0.004, 0.004, -0.000697283865190439, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, -0.000915431872142883, 0.00923966021143272, 0.0056517320473431, 0.00433388821481831, 0.00318723142202518, 0.00340970571823595, 0.00422203922755135, 0.00453976680412294, 0.0046202877861286, 0.00473675955977054, -0.00015508535585088, 0.0101430947793994, -0.000200687523000682, 0.004, 0.004, 0.00192670843135786, -0.000170668681166593, 0.004, 0.004, 0.004, 0.00134722059955132, -0.000262884611492943, 0.004, 0.004, 0.004, 0.00591055345048618, -0.00194419654256857, 0.004, 0.004, 0.004, 0.0006059521491806, -0.000666125818866443, 0.004, 0.004, 0.004, -0.000848820207580652, -0.00202603929390789, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, -0.00421963564488254, 0.0197061148874491, -0.000264471542406444, 0.00285455039246963, 0.00430789439341713, 0.00497175435171098, 0.00520666482053443, 0.00512970956948983, 0.00756654678975618, -0.000156210177176093, -0.000167385712828239, 0.0235174249487206, -0.00446261437531541, 0.004, -0.000137292454515056, -0.000113727761598494, 0.004, 0.004, 0.004, 0.00195878853690182, -0.00243036614907955, 0.004, 0.004, 0.004, -0.00012589615196831, -0.00149689910704691, 0.004, 0.004, 0.004, 0.00292208772529403, -0.000992719687024317, 0.004, 0.004, 0.004, -0.000578660129118478, -0.00154971276165391, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.00491898584259906, 0.00276870412397827, 0.00494455774812606, 0.00656766760623052, 0.00727402059485561, 0.00790510525317462, 0.00695734080343227, 0.00795911616244319, 0.00149771814724309, 0.0064872647513368, 0.00838941933783797, 0.0503795079249428, 0.133565403436032, 0.00118111799840934, -0.000111966052028284, 0.004, 0.004, 0.004, 0.00297489888596317, -0.000616551518766928, 0.004, 0.004, 0.004, 0.00469478817973125, -0.000769862434642523, 0.004, 0.004, 0.004, -0.00015699987838295, 0.0088374935085937, 0.004, 0.004, 0.004, 0.0202955256196721, -0.00152105472223549, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, -0.00162485948842415, 0.00719031288216238, 0.0117184425418287, 0.00650778161041364, 0.00520731821545655, 0.00626210936476919, 0.00528387320775105, 0.00241329476750947, 0.00222541926230401, 0.00486105177918752, -0.000137992232836665, 0.0551088997284223, 0.00277581090829659, -0.000130489857128356, 0.004, 0.004, 0.004, -0.00105027521524037, 0.00421240391851075, 0.004, 0.004, 0.004, -0.000246495319370642, 0.00752914495257707, 0.004, 0.004, 0.004, -0.000455549647551561, 0.00391330861938618, 0.004, 0.004, 0.004, -0.00506357089028416, -0.00112587601697093, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.394658127207145, 0.132198623338349, 0.122779672297155, -0.000787701729946966, -0.000136773490015179, -0.000326350252038345, -0.000160577124974238, 0.00588741182298707, -3.31833880405605E-05, 0.0659420306907641, 0.119489865608271, 0.000906659443464936, -0.000460157995175247, 0.004, 0.004, 0.004, 0.00453739023451891, -0.00100206634666101, 0.004, 0.004, 0.004, 0.00403250216799566, -0.00119057735244476, 0.004, 0.004, 0.004, 0.00277244308388306, -0.000127162097771313, 0.004, 0.004, 0.004, -0.000168637450522381, -0.000698080659819837, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, -0.000379259851061545, 0.008930926532953, 0.00705386419964922, 0.00368061212439558, 0.00349628921563389, 0.00184085339117002, 0.00252589965226819, 0.00297244615576368, 0.0148402531245904, 0.0175789399147927, 0.022328701157952, 0.0188790551362496, -0.000253005586834373, 0.004, 0.004, 0.000842464090188183, -0.000417509287812737, 0.004, 0.004, 0.004, -0.000219265001748337, -0.000642050739424298, 0.004, 0.004, 0.004, -0.000240190513843333, -0.00243004105821436, 0.004, 0.004, 0.004, -0.000338893840987368, -0.00391503152590853, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, -0.109993601332107, -0.109993601311115, -0.109993601459089, -0.0133670221877993, 0.00232441190071383, -0.000128105395900993, 0.0107338200028579, -0.000130417550122584, 0.138246889117962, -0.00013023903862643, 0.020993427908333, 0.0230762421048228, 0.172880843964545, 0.004, -0.00348352964218206, 0.00411608259520258, 0.004, 0.004, 0.004, 0.0240751417424512, -0.000487442887420718, 0.004, 0.004, 0.004, -0.000607923747538089, -0.000740448669514782, 0.004, 0.004, 0.004, -0.000643336986985209, -0.000158910769727784, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 2.1227776026512, 2.25992582522148, 2.06331304123723, -0.000199345727005251, -0.000135482960452555, -0.00381441967855232, -0.000157566209191, -0.000952614525146895, 0.00789035617904118, -0.00029433404969332, 0.0121285401742267, 0.138783444344786, 0.0529849519358424, 0.00860490424520297, 0.00658248537022571, 0.004, 0.004, 0.004, 0.0181498454883923, -0.000186243234500176, 0.004, 0.004, 0.004, 0.0161619894432673, -0.000116094266358386, 0.004, 0.004, 0.004, 0.263742291098308, -0.00254706434660157, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.0204331793450363, 0.0240088328087552, 0.0118546874297609, 0.00645316827422745, 0.0300920966583917, 0.100868412150879, 0.230459456758348, 0.034301607706075, 0.00787222228670539, 0.0322862931710811, 0.107292222613912, 0.103550562241627, -0.000132681071538273, -0.000922380548960452, 0.004, 0.004, 0.004, -0.00287874326257359, 0.00099663634834937, 0.004, 0.004, 0.004, -0.00128150821598214, 0.00570938044135608, 0.004, 0.004, 0.004, -0.00208780047793811, -0.00180873055592343, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, -0.000239936641638167, -0.000287603387052978, -0.00429263576586227, -0.000371451495500384, -0.00184959318433455, -0.000142838210890787, -0.000179925409573789, 4.68893898482993E-05, -0.000258226041305802, 0.0204605293928246, 0.00618641280580659, 0.0100777004107312, 0.00514965369473394, -0.000295356616052253, 0.00781644091887058, 0.00139384322772456, 0.00577230392273254, 0.00139463045314927, -0.00379759219133755, 0.00311785033020844, -0.000108243724912957, 0.00420175495115992, -0.00016454429682426, -0.00156292493154358, -0.000279465662044163, -0.000107542649539694, 0.0495656136591202, -0.00718763061288769, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.0957602894247488, 0.088505252712156, 0.084539296232394, -0.000326593022090893, 0.00910414712742187, 0.0238844279293658, -0.000214358210562323, 0.0125287260866412, 0.0160594427311772, 0.0117496109222664, 0.0230134245178559, -0.00044087903132657, 0.0462713370533124, -0.000115491702336862, -0.000157492029160331, 0.00308867046785184, -0.000334734023139097, 0.000305104920306903, -0.00442061855034225, -0.000586404099713766, -0.000246446717356828, -0.00274095313668139, -0.000334325688944815, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.0220838810096625, -0.000505982726252018, -0.00077412283945304, -0.00047165995860855, -0.000161014859172399, -0.00545669272342152, 0.0164948802518483, 0.0124610574303309, 0.0303527412281259, 0.046304097235554, -0.000112019943446957, 0.0178629056426194, 0.00779302599967594, 0.0108021572160552, 0.0398937599329156, -0.000215427682184859, 0.0303404569581742, 0.0565081117677152, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.802517658984342, -0.00794287927733691, -0.0442200662873663, -0.000112885343174783, -0.000199713996031243, -0.000211692984338995, 0.000337909961358949, -0.000123876763619046, 0.00134278687038659, -0.000220368768960302, 0.00121832676286868, -0.000138440355803643, 0.0517570250183169, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004 };	
//			//LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTimeHomogenousPiecewiseConstant(timeDiscretizationFromArray, liborPeriodDiscretization, new TimeDiscretizationFromArray(0.00, 0.50, 1.00, 2.00, 3.00, 4.00, 5.00, 7.00, 10.00, 15.00, 20.00, 25.00, 30.00 ), new double[] { 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0, 0.20/100.0 });
//			//LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstant(timeDiscretizationFromArray, liborPeriodDiscretization, new TimeDiscretizationFromArray(0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 11.00, 12.00, 14.00, 15.00, 16.00, 17.00, 18.50, 20.00, 22.50, 25.00, 30.00 ), new TimeDiscretizationFromArray(0.00, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 11.00, 12.00, 14.00, 15.00, 16.00, 17.00, 18.50, 20.00, 22.50, 25.00, 30.00 ), 0.40 / 100);
//			final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstant(timeDiscretizationFromArray, liborPeriodDiscretization, optionMaturityDiscretization, timeToMaturityDiscretization, 0.40 / 100, true);
		
		 // If simulation time is below libor time, exceptions will be hard to track.
			final double lastTime	= 21.0;
			final double dt		= 0.005;
			final TimeDiscretizationFromArray timeDiscretizationFromArray = new TimeDiscretizationFromArray(0.0, (int) (lastTime / dt), dt);
			final TimeDiscretization liborPeriodDiscretization = timeDiscretizationFromArray;
			
			final BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(timeDiscretizationFromArray, numberOfFactors, numberOfPaths, 31415 /* seed */);
	
			TimeDiscretization timeToMaturityDiscretization = new TimeDiscretizationFromArray(0.00, 1.00, 2.00, 3.00, 4.00, 5.00,  7.00,  10.0, 15.0, 21.0);		// needed if you use LIBORVolatilityModelPiecewiseConstantWithMercurioModification: TimeDiscretization  = new TimeDiscretizationFromArray(0.0, 0.25, 0.50, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0);

			//TimeDiscretization timeToMaturityDiscretization = new TimeDiscretizationFromArray(0.00,	1,0, 2.00, 3.00, 20.0);
			double[] arrayValues = new double [timeToMaturityDiscretization.getNumberOfTimes()];
			for (int i=0; i<timeToMaturityDiscretization.getNumberOfTimes(); i++) {arrayValues[i]= 0.2/100;}
			
			//final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTimeHomogenousPiecewiseConstant(timeDiscretizationFromArray, liborPeriodDiscretization, timeToMaturityDiscretization, arrayValues);
			//final LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, 0.05, false);
			double c=0.3;
			LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialForm(timeDiscretizationFromArray, liborPeriodDiscretization, -0.3/100.0, 0.2/100.0, c, 0.6/100.0 , true); //-0.3/100.0, 0.2/100.0, c, 0.6/100.0   

			System.out.println("speed test2:" + c );
			final LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, 0.05, false);

			//AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelExponentialForm5Param(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, new double[] { 0.20/100.0, 0.05/100.0, 0.10, 0.05/100.0, 0.10} );
			AbstractLIBORCovarianceModelParametric	covarianceModelParametric = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretizationFromArray, liborPeriodDiscretization, volatilityModel, correlationModel);
//
//			final AbstractLIBORCovarianceModelParametric covarianceModelBlended = new BlendedLocalVolatilityModel(covarianceModelParametric, 0.0, false);
//			final AbstractLIBORCovarianceModelParametric covarianceModelDisplaced = new DisplacedLocalVolatilityModel(covarianceModelParametric, 1.0/0.25, false /* isCalibrateable */);
//
//			final TimeDiscretization tenorTimeScalingDiscretization = new TimeDiscretizationFromArray(0.0, 5.0, 0.5, ShortPeriodLocation.SHORT_PERIOD_AT_START);
//			final double[] tenorTimeScalings = new double[tenorTimeScalingDiscretization.getNumberOfTimes()];
//			Arrays.fill(tenorTimeScalings, 0.0);
//			final TermStructureTenorTimeScalingInterface tenorTimeScalingModel = new TermStructureTenorTimeScalingPicewiseConstant(tenorTimeScalingDiscretization, tenorTimeScalings);

			TermStructureCovarianceModelParametric  termStructureCovarianceModel = new TermStructCovarianceModelFromLIBORCovarianceModelParametric (null, covarianceModelParametric );
//			termStructureCovarianceModel = termStructureCovarianceModel.getCloneWithModifiedParameters(scalingAndVol);

//				System.out.println("Number of volatility parameters: " + volatilityModel.getParameter().length);
//				System.out.println("Number of scaling parameters: " + tenorTimeScalingModel.getParameter().length);
				System.out.println("Number of covariance + scaling parameters: " + termStructureCovarianceModel.getParameter().length);
			
				/*
				 * Create corresponding LIBOR Market Model
				 * 7 days + 3 weeks + 2 month, 1 quarter + 39 semi-annuals
				 */
				final CalibrationProduct[] calibrationItems = new CalibrationProduct[0];

				final TimeDiscretization liborPeriodDiscretizationDaily = new TimeDiscretizationFromArray(0.0, 21.0, 0.005, ShortPeriodLocation.SHORT_PERIOD_AT_START);
				final TimeDiscretization liborPeriodDiscretizationWeekly = new TimeDiscretizationFromArray(0.0, 21.0, 0.025, ShortPeriodLocation.SHORT_PERIOD_AT_START);
				final TimeDiscretization liborPeriodDiscretizationMonthly = new TimeDiscretizationFromArray(0.0, 21.0, 0.125, ShortPeriodLocation.SHORT_PERIOD_AT_START);
				final TimeDiscretization liborPeriodDiscretizationQuarterly = new TimeDiscretizationFromArray(0.0, 21.0, 0.25, ShortPeriodLocation.SHORT_PERIOD_AT_START);
				final TimeDiscretization liborPeriodDiscretizationSemiannual = new TimeDiscretizationFromArray(0.0, 21.0, 0.5, ShortPeriodLocation.SHORT_PERIOD_AT_START);

				final TermStructureModel  liborMarketModelNOTCalibrated = new LIBORMarketModelWithTenorRefinement (
						new TimeDiscretization[] { liborPeriodDiscretizationDaily, liborPeriodDiscretizationWeekly, liborPeriodDiscretizationMonthly,liborPeriodDiscretizationQuarterly,liborPeriodDiscretizationSemiannual},
						new Integer[] { 5, 4, 3, 2, 200 },
						curveModel,
						forwardCurve,
						new DiscountCurveFromForwardCurve(forwardCurve),
						termStructureCovarianceModel,
						calibrationItems,
						null		
				);

				final EulerSchemeFromProcessModel process = new EulerSchemeFromProcessModel(liborMarketModelNOTCalibrated,brownianMotion);
				LIBORModelMonteCarloSimulationModel  simulationNONCalibrated = new LIBORMonteCarloSimulationFromTermStructureModel ( process);

		
//		// CAPLET ON BACKWARD LOOKING RATE SEMESTRALI
		DecimalFormat formatterTimeValue = new DecimalFormat("##0.00;");
		DecimalFormat formatterVolValue = new DecimalFormat("##0.00000;");
		DecimalFormat formatterAnalytic = new DecimalFormat("##0.000;");
		DecimalFormat formatterPercentage = new DecimalFormat(" ##0.000%;-##0.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	
		double strike = 0.004783;
		double dtLibor = 0.5;
		
		double[] mktData = new double[] {/* 6M 0.00167, */ /* 12M*/ 0.00201, /* 18M*/ 0.00228, /* 2Y */ 0.00264, 0.0, /* 3Y */ 0.0033, /* 4Y */0.00406, /* 5Y */ 0.00455, /* 6Y - NA */ 0.0, /* 7Y */0.00513, /* 8Y- NA */0.0, /* 9Y */0.0, /* 10Y */0.00550,0.0,0.0,0.0,0.0, /* 15Y */0.00544,0.0,0.0,0.0,0.0, /* 20Y */0.0053,0.0,0.0,0.0};
	
		int mktDataIndex = 0;
		double beginLiborTime = 0.5;

		System.out.println("\n results on NON CALIBRATED model \n");
		

		while(beginLiborTime < liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1)) {
			Caplet capletCassical = new Caplet(beginLiborTime, dtLibor, strike, dtLibor, false, ValueUnit.NORMALVOLATILITY);
			CapletOnCompoundedBackward backwardCaplet = new CapletOnCompoundedBackward(beginLiborTime, dtLibor, strike, dtLibor,false);
			double capletMaturity = beginLiborTime+dtLibor;
			double impliedVolClassical = capletCassical.getValue(simulationNONCalibrated);
			double impliedVolBackward = backwardCaplet.getValue(simulationNONCalibrated);
			double analyticFormulaPaper = Math.sqrt(1+0.5/(beginLiborTime*3));
			double ratioImpliedVol = impliedVolBackward/impliedVolClassical;
			double error = (analyticFormulaPaper-ratioImpliedVol)/analyticFormulaPaper;

			if (beginLiborTime<2.5) {		//da i valori del caplet per maturity 1.0, 1.5Y, 2Y, 2.5Y, 3.0Y.
				if (mktData[mktDataIndex] == 0.0) {
				System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
						+ formatterTimeValue.format(capletMaturity) + "; "
						+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
						+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
						+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
						+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
						+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
					+ formatterPercentage.format(error) );
				beginLiborTime+=0.5;
				mktDataIndex+=1;
				}
				else {

					double ratioMktVol =impliedVolBackward/mktData[mktDataIndex];
					double errorMkt = (analyticFormulaPaper-ratioMktVol)/analyticFormulaPaper;
					System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
							+ formatterTimeValue.format(capletMaturity) + "; "
							+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
							+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
							+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
							+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
							+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
							+ formatterPercentage.format(error) 
							+  "\t" +"Market Caplet Vol. " + formatterPercentage.format(mktData[mktDataIndex])
							+  "\t" +"Market Ratio: "	+ formatterAnalytic.format(ratioMktVol)	
							+  "\t" +"Market Error: "	+ formatterPercentage.format(errorMkt)
							);
					beginLiborTime+=0.5;
					mktDataIndex+=1;
				}
				
			}
			else {	//secondo loop da i valori del caplet per maturity 4Y,5Y,...,21
	
				if (mktData[mktDataIndex] == 0.0) {
					System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
							+ formatterTimeValue.format(capletMaturity) + "; "
							+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
							+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
							+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
							+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
							+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
						+ formatterPercentage.format(error) );
					beginLiborTime+=1;
					mktDataIndex+=1;
					}
					else {
						double ratioMktVol =impliedVolBackward/mktData[mktDataIndex];
						double errorMkt = (analyticFormulaPaper-ratioMktVol)/analyticFormulaPaper;
						System.out.println("Caplet on B(" + formatterTimeValue.format(beginLiborTime) + ", "
								+ formatterTimeValue.format(capletMaturity) + "; "
								+ formatterTimeValue.format(capletMaturity) + ")." + "\t" + "Backward caplet: "
								+ formatterPercentage.format(impliedVolBackward)+ "\t" + "Classical Caplet: " 
								+ formatterPercentage.format(impliedVolClassical)+ "\t" + "Analytic ratio: " 
								+ formatterAnalytic.format(analyticFormulaPaper) +  "\t" + "MC ratio: "
								+ formatterAnalytic.format(ratioImpliedVol) +  "\t" +"Error: "
								+ formatterPercentage.format(error) 
								+  "\t" +"Market Caplet Vol. " + formatterPercentage.format(mktData[mktDataIndex])
								+  "\t" +"Market Ratio: "	+ formatterAnalytic.format(ratioMktVol)	
								+  "\t" +"Market Error: "	+ formatterPercentage.format(errorMkt)
								);
						beginLiborTime+=1;
						mktDataIndex+=1;
					}
			}
		}
		
		System.out.println("\n Backward looking rate:");
		for(beginLiborTime=0.0; beginLiborTime<liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1); beginLiborTime+=1) {
			double endLiborTime=beginLiborTime+dtLibor;
			RandomVariable backwardLookingRate =  simulationNONCalibrated.getBackward(endLiborTime, beginLiborTime,endLiborTime);
			double avgBackwardLookingRate =backwardLookingRate.getAverage();
			System.out.println("Backward B(" + formatterTimeValue.format(beginLiborTime) +  ", " + formatterTimeValue.format(endLiborTime) + "; "  + formatterTimeValue.format(endLiborTime)+ ")."+ "\t" + "Average: " + avgBackwardLookingRate);
		}
		
		System.out.println("\n LIBOR rate:");
		for(beginLiborTime=0.0; beginLiborTime<liborPeriodDiscretization.getTime(liborPeriodDiscretization.getNumberOfTimes()-1); beginLiborTime+=1) {
			double endLiborTime=beginLiborTime+dtLibor;
			RandomVariable libor =  simulationNONCalibrated.getLIBOR(beginLiborTime, beginLiborTime,endLiborTime);
			double avgBackwardLookingRate = libor.getAverage();
			System.out.println("Backward B(" + formatterTimeValue.format(beginLiborTime) +  ", " + formatterTimeValue.format(endLiborTime) + "; "  + formatterTimeValue.format(beginLiborTime)+ ")."+ "\t" + "Average: " + avgBackwardLookingRate);
		}
		
}
		
		public AnalyticModel getCalibratedCurve() throws SolverException {
			final String[] maturity					= { "1D", "7D", "14D", "21D", "1M", "2M", "3M", "4M", "5M", "6M", "7M", "8M", "9M", "12M", "15M", "18M", "21M", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "11Y", "12Y", "15Y", "20Y", "25Y", "30Y", "40Y", "50Y" };
			final String[] frequency				= { "tenor", "tenor", "tenor",  "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual"};
			final String[] frequencyFloat			= { "tenor", "tenor", "tenor",  "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual"};
			final String[] daycountConventions	    = { "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360"};
			final String[] daycountConventionsFloat	= { "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360"};
			final double[] rates					= { -0.0055, -0.00553, -0.00553, -0.00553, -0.00553, -0.00555, -0.00556, -0.00559, -0.00564, -0.00568, -0.00572, -0.00577, -0.00581, -0.00592, -0.00601, -0.00608, -0.00613, -0.00619, -0.00627, -0.00622, -0.00606, -0.00582, -0.00553, -0.00519, -0.00482, -0.00442, -0.00402, -0.00362, -0.00261, -0.00189, -0.00197, -0.0023, -0.00286, -0.00333};
			final HashMap<String, Object> parameters = new HashMap<>();

			parameters.put("referenceDate", LocalDate.of(2020, Month.JULY, 31));
			parameters.put("currency", "EUR");
			parameters.put("forwardCurveTenor", "1D");
			parameters.put("maturities", maturity);
			parameters.put("fixLegFrequencies", frequency);
			parameters.put("floatLegFrequencies", frequencyFloat);
			parameters.put("fixLegDaycountConventions", daycountConventions);
			parameters.put("floatLegDaycountConventions", daycountConventionsFloat);
			parameters.put("rates", rates);

		return getCalibratedCurve(null, parameters);
	}

	private static AnalyticModel getCalibratedCurve(final AnalyticModel model2, final Map<String, Object> parameters) throws SolverException {

		final LocalDate	referenceDate		= (LocalDate) parameters.get("referenceDate");
		final String	currency			= (String) parameters.get("currency");
		final String	forwardCurveTenor	= (String) parameters.get("forwardCurveTenor");
		final String[]	maturities			= (String[]) parameters.get("maturities");
		final String[]	frequency			= (String[]) parameters.get("fixLegFrequencies");
		final String[]	frequencyFloat		= (String[]) parameters.get("floatLegFrequencies");
		final String[]	daycountConventions	= (String[]) parameters.get("fixLegDaycountConventions");
		final String[]	daycountConventionsFloat	= (String[]) parameters.get("floatLegDaycountConventions");
		final double[]	rates						= (double[]) parameters.get("rates");

		Assert.assertEquals(maturities.length, frequency.length);
		Assert.assertEquals(maturities.length, daycountConventions.length);
		Assert.assertEquals(maturities.length, rates.length);

		Assert.assertEquals(frequency.length, frequencyFloat.length);
		Assert.assertEquals(daycountConventions.length, daycountConventionsFloat.length);

		final int		spotOffsetDays = 2;
		final String	forwardStartPeriod = "0D";

		final String curveNameDiscount = "discountCurve-" + currency;

		/*
		 * We create a forward curve by referencing the same discount curve, since
		 * this is a single curve setup.
		 *
		 * Note that using an independent NSS forward curve with its own NSS parameters
		 * would result in a problem where both, the forward curve and the discount curve
		 * have free parameters.
		 */
		final ForwardCurve forwardCurve		= new ForwardCurveFromDiscountCurve(curveNameDiscount, referenceDate, forwardCurveTenor);

		// Create a collection of objective functions (calibration products)
		final Vector<AnalyticProduct> calibrationProducts = new Vector<>();
		final double[] curveMaturities	= new double[rates.length+1];
		final double[] curveValue			= new double[rates.length+1];
		final boolean[] curveIsParameter	= new boolean[rates.length+1];
		curveMaturities[0] = 0.0;
		curveValue[0] = 1.0;
		curveIsParameter[0] = false;
		for(int i=0; i<rates.length; i++) {
//			scheduleRec=scheduleReceiveLeg; schedulePay=schedulePayLeg
			final Schedule schedulePay = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturities[i], frequency[i], daycountConventions[i], "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);
			final Schedule scheduleRec = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturities[i], frequencyFloat[i], daycountConventionsFloat[i], "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);

			curveMaturities[i+1] = Math.max(schedulePay.getPayment(schedulePay.getNumberOfPeriods()-1),scheduleRec.getPayment(scheduleRec.getNumberOfPeriods()-1));
			curveValue[i+1] = 1.0;
			curveIsParameter[i+1] = true;
			calibrationProducts.add(new Swap(schedulePay, null, rates[i], curveNameDiscount, scheduleRec, forwardCurve.getName(), 0.0, curveNameDiscount));
		} 

		final InterpolationMethod interpolationMethod = InterpolationMethod.LINEAR;

		// Create a discount curve
		final DiscountCurveInterpolation discountCurveInterpolation = DiscountCurveInterpolation.createDiscountCurveFromDiscountFactors(
				curveNameDiscount								/* name */,
				referenceDate	/* referenceDate */,
				curveMaturities	/* maturities */,
				curveValue		/* discount factors */,
				curveIsParameter,
				interpolationMethod ,
				ExtrapolationMethod.CONSTANT,
				InterpolationEntity.LOG_OF_VALUE
				);

		/*
		 * Model consists of the two curves, but only one of them provides free parameters.
		 */
		AnalyticModel model = new AnalyticModelFromCurvesAndVols(new Curve[] { discountCurveInterpolation, forwardCurve });

		/*
		 * Create a collection of curves to calibrate
		 */
		final Set<ParameterObject> curvesToCalibrate = new HashSet<>();
		curvesToCalibrate.add(discountCurveInterpolation);

		/*
		 * Calibrate the curve
		 */
		final Solver solver = new Solver(model, calibrationProducts);
		final AnalyticModel calibratedModel = solver.getCalibratedModel(curvesToCalibrate);
		System.out.println("Solver reported acccurary....: " + solver.getAccuracy());

		Assert.assertEquals("Calibration accurarcy", 0.0, solver.getAccuracy(), 1E-3);

		// Get best parameters
		final double[] parametersBest = calibratedModel.getDiscountCurve(discountCurveInterpolation.getName()).getParameter();

		// Test calibration
		model			= calibratedModel;

//--> praticamente il calibrato valore di calibrationProducts (che sono swap) deve essere = 0 --> penso perchè il prezzo dello swap deve essere = 0 per definizione, quindi in teoria gli interest definiti sono L_i - S e quando calibri questa differenza deve essere =0 (e se così è allora il prezzo dello swap è 0) --> avrebbe anche senso perchè se noti il payer ha tassi prima negativi e poi positivi, mentre quell'altro ha tassi sempre = 0 (cioè L=S, cioè è quello che paga il fixed swap rate)
		double squaredErrorSum = 0.0;
		for(final AnalyticProduct c : calibrationProducts) {
			final double value = c.getValue(0.0, model);
			final double valueTaget = 0.0;
			final double error = value - valueTaget;
			squaredErrorSum += error*error;
		}
		final double rms = Math.sqrt(squaredErrorSum/calibrationProducts.size());

		System.out.println("Independent checked acccurary: " + rms);

		System.out.println("Calibrated discount curve: ");
		for(int i=0; i<curveMaturities.length; i++) {
			final double maturity = curveMaturities[i];
			System.out.println(maturity + "\t" + calibratedModel.getDiscountCurve(discountCurveInterpolation.getName()).getDiscountFactor(maturity));
		}
		return model;
	}

	private static double getParSwaprate(final ForwardCurve forwardCurve, final DiscountCurve discountCurve, final double[] swapTenor) {
		return net.finmath.marketdata.products.Swap.getForwardSwapRate(new TimeDiscretizationFromArray(swapTenor), new TimeDiscretizationFromArray(swapTenor), forwardCurve, discountCurve);
	}
}
