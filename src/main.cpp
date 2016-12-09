/*-
 * The MIT License (MIT)
 *
 * Copyright (c) 2014, 2015 Guram Duka
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
//------------------------------------------------------------------------------
#include "nn.hpp"
//------------------------------------------------------------------------------
int main()
{
	using namespace std;
	using namespace nn;

	cout.precision(120);

#if _WIN32
	//SetPriorityClass(GetCurrentProcess(), REALTIME_PRIORITY_CLASS);
	//SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
	//SetProcessAffinityMask(GetCurrentProcess(), 1);
	//SetThreadAffinityMask(GetCurrentThread(), 1);

	LARGE_INTEGER freq, start, stop, hole, hole_start;
	QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&start);
	hole.QuadPart = 0;
#else
    struct timespec res, start, stop, hole_start;
    const uint64_t freq = 1000000000u;
    uint64_t hole = 0;

    auto ts2i = [] (struct timespec & ts) {
        return uint64_t(ts.tv_sec) * freq + ts.tv_nsec;
    };

    clock_getres(CLOCK_PROCESS_CPUTIME_ID,&res);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start);
#endif

	numeric hard, soft, error;
#if 0
	{
		random_device rd;
		default_random_engine e1(rd());
		uniform_int_distribution<int> uniform_dist('0', '9');
		string s;

		for( int i = 0; i < 10000; i++ )
			s.push_back(uniform_dist(e1));
		s.push_back('.');
		for( int i = 0; i < 10000; i++ )
			s.push_back(uniform_dist(e1));
		hard = numeric(s);
		hard.normalize(3);
		QueryPerformanceCounter(&start);
		s = hard.to_string(0, 10000);
		LARGE_INTEGER stop;
		QueryPerformanceCounter(&stop);

		uint64_t ellapsed  = stop.QuadPart - start.QuadPart;
		uint64_t milisecs  = (ellapsed * 1000u) / freq.QuadPart;
		uint64_t microsecs = (ellapsed * 1000000u) / freq.QuadPart;
		uint64_t nanosecs  = (ellapsed * 1000000000u) / freq.QuadPart;
		cout << "MILIS  : " << milisecs << endl;
		cout << "MICROS : " << microsecs << endl;
		cout << "NANOS  : " << nanosecs << endl;

		stringstream ss;
		ss
			<< s << ", stat: " << hard.stat_length() << endl;
	}
	//return 0;
#endif

	hard = numeric(
		// http://oeis.org/A002193
		// http://apod.nasa.gov/htmltest/gifcity/sqrt2.1mil
		"1.4142135623730950488016887242096980785696718753769480731766797379907324784621070"
		"388503875343276415727350138462309122970249248360558507372126441214970999358314132"
		"226659275055927557999505011527820605714701095599716059702745345968620147285174186"
		"408891986095523292304843087143214508397626036279952514079896872533965463318088296"
		"406206152583523950547457502877599617298355752203375318570113543746034084988471603"
		"868999706990048150305440277903164542478230684929369186215805784631115966687130130"
		"156185689872372352885092648612494977154218334204285686060146824720771435854874155"
		"657069677653720226485447015858801620758474922657226002085584466521458398893944370"
		"926591800311388246468157082630100594858704003186480342194897278290641045072636881"
		"313739855256117322040245091227700226941127573627280495738108967504018369868368450"
		"725799364729060762996941380475654823728997180326802474420629269124859052181004459"
		"842150591120249441341728531478105803603371077309182869314710171111683916581726889"
		"419758716582152128229518488472");
	hard.normalize(3);
	soft = nn::numeric(2).root(2, 10);
	soft.normalize(3);
	error = (soft - hard).abs();
	error.normalize(3);

#if _WIN32
	QueryPerformanceCounter(&hole_start);
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&hole_start);
#endif
	cout
		<< "sqrt(2):" << endl
		<< "hard   : " << hard.to_string(0, 1000) << ", stat: " << hard.stat_length() << endl
		<< "soft   : " << soft.to_string(0, 1000) << ", stat: " << soft.stat_length() << endl
		<< "error  : " << error.to_string(0, 1000) << ", stat: " << error.stat_length() << endl;
#if _WIN32
	QueryPerformanceCounter(&stop);
	hole.QuadPart += stop.QuadPart - hole_start.QuadPart;
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop);
    hole += ts2i(stop) - ts2i(hole_start);
#endif

	hard = numeric(
		"1.2599210498948731647672106072782283505702514647015079800819751121552996765139594"
		"837293965624362550941543102560356156652593990240406137372284591103042693552469606"
		"426166250009774745265654803068671854055186892458725167641993737096950983827831613"
		"991551293136953661839474634485765703031190958959847411059811629070535908164780114"
		"735213254847712978802422085820532579725266622026690056656081994715628176405060664"
		"826773572670419486207621442965694205079319172441480920448232840127470321964282081"
		"201905714188996459998317503801888689594202055922021154729973848802607363697417887"
		"792157984675099539630078260959624203483238660139857363433909737126527995991969968"
		"377913168168154428850279651529278107679714002040605674803938561251718357006907984"
		"996341976291474044834540269715476228513178020643878047649322579052898467085805286"
		"258130005429388560720609747223040631357234936458406575916916916727060124402896700"
		"001069081035313852902700415084");
	hard.normalize(3);
	soft = nn::numeric(2).root(3, 5100);
	soft.normalize(3);
	error = (soft - hard).abs();

#if _WIN32
	QueryPerformanceCounter(&hole_start);
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&hole_start);
#endif
	cout
		<< "cube_root(2):" << endl
		<< "hard   : " << hard.to_string(0, 1000) << ", stat: " << hard.stat_length() << endl
		<< "soft   : " << soft.to_string(0, 1000) << ", stat: " << soft.stat_length() << endl
		<< "error  : " << error.to_string(0, 1000) << ", stat: " << error.stat_length() << endl;
#if _WIN32
	QueryPerformanceCounter(&stop);
	hole.QuadPart += stop.QuadPart - hole_start.QuadPart;
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop);
    hole += ts2i(stop) - ts2i(hole_start);
#endif

	hard = numeric(
		// computed by http://pythonhosted.org/bigfloat/
		"0.841470984807896506652502321630298999622563060798371065672751709991910404391239"
		"66894863974354305269585434903790792067429325911892099189888119341032772921240948"
		"07919558267666069999077640119784087827325663474848028702986561570179624553948935"
		"72924670127086486281053382030561377218203868449667761674266239013382753397956764"
		"25556547796398976482432869027569642912063005830365152303127825528985326485139819"
		"34521359709559620621721148144417810576010756741366480550089167266058041400780623"
		"93070371877956261288804636081734524656391420252404187763420749206952007713347809"
		"81427902145268255663208233521544160916442090589298702247338446044897237139799127"
		"40819247250488554873119310350681908151532607457392911183319628215089734868811421"
		"45283822986512570166738407445519237561432212906059248273970368180158563090543266"
		"78464310753126381217325670198560110683602890189501942151616655191791451720046686"
		"595971691072197805885406460019");
	hard.normalize(3);
	soft = nn::numeric(1).sin(224);
	soft.normalize(3);
	error = (soft - hard).abs();

#if _WIN32
	QueryPerformanceCounter(&hole_start);
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&hole_start);
#endif
	cout
		<< "sin(1) :" << endl
		<< "hard   : " << hard.to_string(0, 1000) << ", stat: " << hard.stat_length() << endl
		<< "soft   : " << soft.to_string(0, 1000) << ", stat: " << soft.stat_length() << endl
		<< "error  : " << error.to_string(0, 1000) << endl;
#if _WIN32
	QueryPerformanceCounter(&stop);
	hole.QuadPart += stop.QuadPart - hole_start.QuadPart;
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop);
    hole += ts2i(stop) - ts2i(hole_start);
#endif

	hard = numeric(
		"0.5403023058681397174009366074429766037323104206179222276700972553811003947744717"
		"645179518560871830893435717311600300890978606337600216634564065122654173185847179"
		"711644744794942331179245513932543359435177567028925963757361543275496417544917751"
		"151312227301006313570782322367714015174689959366787306742276202450776374406758749"
		"816178427202164558511156329688905710812427293316986852471456894904342375433094423"
		"024093596239583182454728173664078071243433621748100322027129757882291764468359872"
		"699426491344391826569453515750762782513804991607306380317214450349861294883363356"
		"557799097930152879278840389800974548251049924537987740061453776371387833594234524"
		"168164283618828482374896327390556260912017589827502528599917438580692485584232217"
		"826858271088291564683006796875955130036108120336747472749181033673515093458888304"
		"203217596594052703934762502487370752661313369842416059710595606599978691384415574"
		"414466420012839398870926323453338868626299654709768054836830358211823411732418465"
		"771864116514294188326444690783");
	hard.normalize(3);
	soft = nn::numeric(1).cos(224);
	soft.normalize(3);
	error = (soft - hard).abs();

#if _WIN32
	QueryPerformanceCounter(&hole_start);
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&hole_start);
#endif
	cout
		<< "cos(1) :" << endl
		<< "hard   : " << hard.to_string(0, 1000) << ", stat: " << hard.stat_length() << endl
		<< "soft   : " << soft.to_string(0, 1000) << ", stat: " << soft.stat_length() << endl
		<< "error  : " << error.to_string(0, 1000) << ", stat: " << error.stat_length() << endl;
#if _WIN32
	QueryPerformanceCounter(&stop);
	hole.QuadPart += stop.QuadPart - hole_start.QuadPart;
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop);
    hole += ts2i(stop) - ts2i(hole_start);
#endif

	// http://oeis.org/A000796
	// http://web.archive.org/web/20140225153300/http://www.exploratorium.edu/pi/pi_archive/Pi10-6.html
	string pi_txt(
		"3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"
		"8214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196"
		"4428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273"
		"7245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094"
		"3305727036575959195309218611738193261179310511854807446237996274956735188575272489122793818301194912"
		"9833673362440656643086021394946395224737190702179860943702770539217176293176752384674818467669405132"
		"0005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235"
		"4201995611212902196086403441815981362977477130996051870721134999999837297804995105973173281609631859"
		"5024459455346908302642522308253344685035261931188171010003137838752886587533208381420617177669147303"
		"5982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"
		"3809525720106548586327886593615338182796823030195203530185296899577362259941389124972177528347913151"
		"5574857242454150695950829533116861727855889075098381754637464939319255060400927701671139009848824012"
		"8583616035637076601047101819429555961989467678374494482553797747268471040475346462080466842590694912"
		"9331367702898915210475216205696602405803815019351125338243003558764024749647326391419927260426992279"
		"6782354781636009341721641219924586315030286182974555706749838505494588586926995690927210797509302955"
		"3211653449872027559602364806654991198818347977535663698074265425278625518184175746728909777727938000"
		"8164706001614524919217321721477235014144197356854816136115735255213347574184946843852332390739414333"
		"4547762416862518983569485562099219222184272550254256887671790494601653466804988627232791786085784383"
		"8279679766814541009538837863609506800642251252051173929848960841284886269456042419652850222106611863"
		"0674427862203919494504712371378696095636437191728746776465757396241389086583264599581339047802759009"
		"9465764078951269468398352595709825822620522489407726719478268482601476990902640136394437455305068203"
		"4962524517493996514314298091906592509372216964615157098583874105978859597729754989301617539284681382"
		"6868386894277415599185592524595395943104997252468084598727364469584865383673622262609912460805124388"
		"4390451244136549762780797715691435997700129616089441694868555848406353422072225828488648158456028506"
		"0168427394522674676788952521385225499546667278239864565961163548862305774564980355936345681743241125");

	numeric pi_v1000(pi_txt), pi_c, pi_e;
	uintptr_t pi_iter = 1, pre = pi_txt.size() - 2;
	pi_cont pi_cnt;

	error = numeric(1) / numeric(10).pow(pre);

	pi_v1000.normalize(3);
	error.normalize(3);

	for(;;){
		pi_c = pi(pi_iter, &pi_cnt);
		pi_e = (pi_v1000 - pi_c).abs();
		if( pi_e < error )
			break;
		if( pi_iter == 1 || pi_iter % 250 == 0 ){
#if _WIN32
			QueryPerformanceCounter(&hole_start);
#else
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&hole_start);
#endif
			cout << "PI  = " << pi_c.to_string(0, pre + 2).c_str()
				<< ", stat: " << pi_c.stat_length();
			cout << ", norm: " << pi_c.normalize(3).stat_length() << endl;
			cout << "err = " << pi_e.to_string(0, pre + 2).c_str() << ", stat: " << pi_e.stat_length() << endl;
#if _WIN32
			QueryPerformanceCounter(&stop);
			hole.QuadPart += stop.QuadPart - hole_start.QuadPart;
#else
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop);
            hole += ts2i(stop) - ts2i(hole_start);
#endif
		}
		pi_iter++;
	}

	pi_c.normalize(3);
	pi_e.normalize(3);

#if _WIN32
	QueryPerformanceCounter(&hole_start);
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&hole_start);
#endif
	cout << "PI(prec) = " << pi_v1000.to_string(0, pre + 2).c_str() << ", stat: " << pi_v1000.stat_length() << endl;
	cout << "PI(calc) = " << pi_c.to_string(0, pre + 2).c_str() << ", stat: " << pi_c.stat_length() << endl;
	cout << "error    = " << pi_e.to_string(0, pre + 2).c_str()  << ", stat: " << pi_e.stat_length() << endl;
	cout << "iterations: " << pi_iter << ", precision: " << pre << endl;
#if _WIN32
	QueryPerformanceCounter(&stop);
	hole.QuadPart += stop.QuadPart - hole_start.QuadPart;
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop);
    hole += ts2i(stop) - ts2i(hole_start);
#endif

#if _WIN32
	QueryPerformanceCounter(&stop);

	uint64_t ellapsed  = stop.QuadPart - start.QuadPart - hole.QuadPart;
	uint64_t milisecs  = (ellapsed * 1000u) / freq.QuadPart;
	uint64_t microsecs = (ellapsed * 1000000u) / freq.QuadPart;
	uint64_t nanosecs  = (ellapsed * 1000000000u) / freq.QuadPart;
#else
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&stop);
    hole += ts2i(stop) - ts2i(hole_start);

	uint64_t ellapsed  = ts2i(stop) - ts2i(start) - hole;
	uint64_t milisecs  = ellapsed * 1000u / freq;
	uint64_t microsecs = ellapsed * 1000000u / freq;
	uint64_t nanosecs  = ellapsed * 1000000000u / freq;
#endif

	cout << endl;
	cout << "SECS   : " << setprecision(2) << std::fixed << (long double) milisecs / 1000u << endl;
	cout << "MILIS  : " << milisecs << endl;
	cout << "MICROS : " << microsecs << endl;
	cout << "NANOS  : " << nanosecs << endl;
	cout << "INDEX  : " <<
		microsecs / (integer::stat_iadd_ + integer::stat_isub_ + integer::stat_imul_ + integer::stat_idiv_)
		<< endl;

	return 0;
}
//------------------------------------------------------------------------------
