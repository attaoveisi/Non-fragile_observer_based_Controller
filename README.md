# Non-fragile_observer_based_Controller
Robust Non-Fragile Observer-based  Controller


A robust non-fragile observer-based controller for linear time-invariant (LTI) system with structured uncertainty is introduced. The  robust stability of the closed-loop system is guaranteed by use of Lyapunov theorem in the presence of undesirable disturbance. For the sake of addressing the fragility problem, independent sets of time dependent gain-uncertainties are assumed to be existing for the controller and the observer elements. In order to satisfy the arbitrary -normed constraints for the control system and to enable automatic determination of the optimal  bound of the performance functions in disturbance rejection control (DRC), additional necessary and sufficient conditions are presented in a linear matrix equality/inequality (LME/LMI) framework. The  observer-based controller is then transformed into an optimization problem of coupled set of LMIs/LME that can be solved iteratively by use of numerical software such as Scilab. Finally, concerning the evaluation of the performance of the controller, the control system is implemented in real-time on a mechanical system with aiming at vibration suppression. The plant under study is a multi-input single-output (MISO) clamped-free piezo-laminated smart beam. The nominal mathematical reduced-order model of the beam with piezo-actuators is used to design the proposed controller and then the control system is implemented experimentally on the full-order real-time system. The results show that the closed-loop system has a robust performance in rejecting the disturbance in the presence of the structured uncertainty and in the presence of the unmodelled dynamics. 

Keywords: Laser velocimetry, Lyapunov methods, Piezoelectric transducers, Robust control, Vibration control.


refer to https://scholar.google.com/citations?user=-HRHoYoAAAAJ&hl=de

References

Dorato P (1998) Non-fragile controller design: and overview. In: Proceedings of the American control conference, Philadelphia, 26-26 June, pp. 2829-2831.

Corrado JR and Haddad WM (1999) Static output feedback controllers for systems with parametric uncertainty and controller gain variation. In: Proceedings of the American control conference, San Diego, California, 2-4 June, pp. 915-919.

Wang R, Xing J, Wang P and Yang Q (2012) Non-fragile observer design for nonlinear switched time delay systems using delta operator. In: International conference control (UKACC), Cardiff, UK, 3-5 September, pp. 387-393. 

Lien C, Cheng W, Tsai C and Yu K (2007) Non-fragile observer-based controls of linear system via LMI approach. Chaos, Solitons & Fractals 32: 1530–1537.

Du H, Lam K and Sze SY (2004) Non-fragile  vibration control for uncertain structural systems. Journal of Sound and Vibration 273: 1031–1045.

Yazici H, Guclu R, Kucukdemiral IB and Parlakci MNA (2012) Robust delay-dependent  control for uncertain structural systems with actuator delay. Journal of Dynamic Systems, Measurement and Control 134: 031013-15.

Ramakrishnan K and Ray G (2012) Delay-dependent non-fragile  robust control for a class of uncertain stochastic nonlinear systems. In: 7th IEEE International Conference on Industrial and Information Systems, Chennai, India.

Takaba K and Katayama T (1998) Robust H2 control of descriptor system with time-varying uncertainty. In: Proceedings of the American control conference, Philadelphia, 26-26 June, pp. 2421-2426.

Yang DM, Zhang QL and Sha CM (2005) Observer-based H2 suboptimal control for descriptor systems. Discrete and Continuous Dynamical Systems, Series B 12: 186-196.

Lai CT, Fang CH, Kau SW and Lee CH (2004) Robust H2 control of norm-bounded uncertain continuous-time system-an LMI approach. In: IEEE International symposium on computer aided control systems design, Taipei, Taiwan, 2-4 September, pp. 243-248.

Yang D, Ma Y, Sha C and Zhang Q (2005) H2 observer for descriptor systems: a LMI design method. International journal of information and systems sciences 1(3-4): 293-301.

Apkarian P, Tuan HD and Bernussou J (2001) Continuous-time analysis, eigenstructure assignment, and H2 synthesis with enhanced linear matrix inequalities (LMI) characterizations. IEEE Transactions on automatic control 46: 1941-1946.

Oveisi A, Gudarzi M, Mohammadi MM and Doosthoseini A (2013) Modeling, identification and active vibration control of a funnel-shaped structure used in MRI throat. Journal of Vibroengineering 15: 1392-8716.

Chen K, Paurobally R, Pan J and Qiu X (2015) Improving active control of fan noise with automatic spectral reshaping for reference signal. Applied acoustics 87: 142-152.

Oveisi A and Gudarzi M (2013) Adaptive sliding mode vibration control of a nonlinear smart beam: A comparison with self-tuning Ziegler-Nichols PID controller. Journal Of Low Frequency Noise Vibration And Active Control 32: 41-62.

Hasheminejad SM and Oveisi A (2015) Active vibration control of an arbitrary thick smart cylindrical panel with optimally placed piezoelectric sensor/actuator pairs. International journal of mechanics and materials in design. DOI: 10.1007/s10999-015-9293-2.

Nestorović T, Durrani N and Trajkov M (2012) Experimental model identification and vibration control of a smart cantilever beam using piezoelectric actuators and sensors. Journal of electroceramics 29: 42-55.

Robu B, Baudouin L, Prieur C and Arzelier D (2012) Simultaneous H∞ vibration control of fluid/plate system via reduced-order controller. IEEE Transactions on control systems technology 20: 3146-3151.

Nestorović T and Trajkov M (2013) Optimal actuator and sensor placement based on balanced reduced models. Mechanical systems and signal processing 36: 271-289.

Oveisi A, Gudarzi M and Hasheminejad SM (2014) Dynamic response of a thick piezoelectric circular cylindrical panel: An exact solution. Shock and vibration, Article ID: 592165.

Boyd S, Ghaoui LE, Feron E and Balakrishan V (1994) Linear matrix inequality in systems and control theory, special for industrial and applied mathematics, SIAM, Philadelphia, PA.

Singh V (2004) Robust stability of cellular neural networks with delay: linear matrix inequality approach. IEE Proceedings-D control theory and applications 151(1): 125–129.

Choi HH (2007) LMI-based sliding surface design for integral sliding mode control of mismatched uncertain systems. IEEE Transactions on automatic control 52: 736-742.

Yang DM, Zhang QL and Yao B (2004) Descriptor Systems, Science Press, Beijing.

Scilab Enterprises, Scilab: Free and Open Source software for numerical computation, 2012.

Soize C (2005) Random matrix theory for modeling uncertainties in computational mechanics. Computer Methods in Applied Mechanics and Engineering, 194(12-16): 1333-66.

Adhikari S, Friswell MI, Lonkar K and Sarkar A (2009) Experimental case studies for uncertainty quantification in structural dynamics. Probabilistic Engineering Mechanics, 24: 473-492.

Takahashi RHC, Dutra DA, Palharess RM and Peresh PLD (2000) On robust non-fragile static state-feedback controller synthesis. In: Proceedings of the 39th IEEE Conference on Decision and Control, Sydney, Australia. 

Duan Z, Huang L and Wang L (2001) Robustness analysis and synthesis of SISO systems under both plant and controller perturbations. Systems & Control Letters 42: 201-216.

Lien C-H, Yu K-W, Lin Y-F, Chung Y-J and Chung L-Y (2008) Robust reliable H∞ control for uncertain nonlinear systems via LMI approach. Applied Mathematics and Computation 198: 453-462. 

Yang G-H and Wang JL (2001) Non-fragile  control for linear systems with multiplicative controller gain variations. Automatica 37: 727-737.

Yang G-H and Wang JL Nonfragile (2003) H∞ output feedback controller design for linear systems. Journal of Dynamic Systems, Measurement, and Control 125(1).

Yang G-H and Che W-W (2006) Non-fragile H∞ filter design with additive gain variations. Proceedings of the 45th IEEE Conference on Decision & Control, Manchester Grand Hyatt Hotel San Diego, CA, USA.

Li G (1998) On the structure of digital controller with finite word length consideration. IEEE Transactions on Automatic Control 43: 689-693.

Liu L, Hana Z and Li W (2010) H∞ non-fragile observer-based sliding mode control for uncertain time-delay systems. Journal of the Franklin Institute 347: 567-576.

Pourgholi M and Majd VJ (2011) A new non-fragile H∞ proportional—integral filtered-error adaptive observer for a class of non-linear systems and its application to synchronous generators. Proceedings of the Institution of Mechanical Engineers, Part I: Journal of Systems and Control Engineering 225(1): 99-112.

Chen J-D, Yang C-D, Lien C-H and Horng J-H (2008) New delay-dependent non-fragile H∞ observer-based control for continuous time-delay systems. Information Sciences 178: 4699-4706. 

Famularo D, Dorato P, Abdallah CT, Haddad MM and Jadbabaie A (2000) Robust non-fragile LQ controllers: The static state feedback case. International Journal of Control 73(2): 159-165.

Changweia Y and Jie C (2012) The Non-fragile Controller Design Based on Quadratic Performance Optimization. AASRI Procedia 3: 2-7. 

Xin-yu M, Oao-bo W and Miao C (2010) Non-fragile Robust Control of System with Coprime Factor Perturbations. International Conference on Computer Application and System Modeling, Taiyuan, china.

Oya H, Hagino K and Mukaidani H (2005) Robust non-fragile controllers for uncertain linear continuous-time systems. 31st Annual Conference of IEEE Industrial Electronics Society, Sheraton Capital Center Raleigh, NC, USA.

Yee J-S, Yang G-H and Wang JL (2001) Non-fragile guaranteed cost control for discrete-time uncertain linear systems. International Journal of Systems Science 32(7): 845-853.

Hsieh C-S, Tyan L, Yang S-S and Lin CC (2005) A resilient guaranteed cost control design via the robust two-stage LQ reliable control. TENCON 2005, Melbourne, Qld, Australia.

Kchaoua M, Souissi M and Toumi A (2012) A new approach to non-fragile H∞ observer-based control for discrete-time fuzzy systems. International Journal of Systems Science 43(1) 9-20.

Tandon A and Dhawan A (2014) An LMI approach to non-fragile robust optimal guaranteed cost control of 2D discrete uncertain systems. Transactions of the Institute of Measurement and Control 36(5): 644-653.

Rajamani R (1998) Observers for Lipschitz Nonlinear Systems. IEEE Transactions on Automatic Control, 43: 397-400.

Mondal S, Chakraborty G and Bhattacharyya K (2010) LMI approach to robust unknown input observer design for continuous systems with noise and uncertainties. International Journal of Control, Automation, and Systems, 8(2): 210-219.

Mahmoodi SN and Jalili N (2007) Non-linear vibrations and frequency response analysis of piezoelectrically driven microcantilevers. International Journal of Non-Linear Mechanics 42(4): 577–587.

Askari H, Esmailzadeh E and Barari A (2015) A unified approach for nonlinear vibration analysis of curved structures using non-uniform rational B-spline representation. Journal of Sound and Vibration 353: 292-307.

Oveisi A and Shakeri R (2016) Robust reliable control in vibration suppression of sandwich circular. Engineering Structures 116: 1–11, DOI: 10.1016/j.engstruct.2016.02.040.

Younesian D, Saadatnia Z and Askari H (2012) Analytical solutions for free oscillations of beams on nonlinear elastic foundations using the variational iteration method. Journal of Theoretical and Applied Mechanics 50(2): 639-652.

Zhang B-L, Meng MM, Han Q-L and Zhang XM (2015) Robust non-fragile sampled-data control for offshore steel jacket platforms. Nonlinear Dynamics 83(4): 1939-1954.

Zhang B-L, Huang ZW and Han QL (2012) Delayed non-fragile H∞ control for offshore steel jacket platforms. In: 38th Annual Conference on IEEE Industrial Electronics Society, Montreal, QC, 25-28 October, pp. 2216-2221.

Zhang B-L, Han Q-L and Huang ZW (2014) Pure delayed non-fragile control for offshore steel jacket platforms subject to nonlinear self-excited wave force. Nonlinear Dynamics 77(3): 491–502.

Oveisi A and Nestorović T (2016) Robust observer-based adaptive fuzzy sliding mode controller. Mechanical Systems and Signal Processing DOI:10.1016/j.ymssp.2016.01.015.
