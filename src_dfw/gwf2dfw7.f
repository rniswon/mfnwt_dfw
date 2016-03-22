      MODULE GWFDFWMODULE
        CHARACTER(LEN=64),PARAMETER:: Version_dfw =
     +'$Id: gwf2dfw7.f 1006 2010-08-15 22:33:45Z Shanafield $'
        DOUBLE PRECISION, SAVE,POINTER:: HEPS                     
        DOUBLE PRECISION, SAVE,POINTER:: CLOSEZERO, QTOL
        INTEGER,SAVE,POINTER:: NSS, ISTCB1, ISTCB2
        INTEGER,SAVE,POINTER:: MAXPTS, IRTFLG
        INTEGER,SAVE,POINTER:: NSTRMdfw, NUMACTIVEDFW
        INTEGER,SAVE,POINTER:: ITMP, IRDFLG, IPTFLG
        REAL   ,SAVE,POINTER:: CONST
        DOUBLE PRECISION,SAVE,POINTER:: TOTSPFLOW
        INTEGER,SAVE,  DIMENSION(:),  POINTER:: IOTSG, NSEGCK,WETTED
        INTEGER,SAVE,  DIMENSION(:,:),POINTER:: ISEGDFW, IDIVA, ISTRMdfw
        REAL,   SAVE,  DIMENSION(:),  POINTER:: SGOTFL
        REAL,   SAVE,  DIMENSION(:),  POINTER:: TSEG,TRCH
        DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER:: QSTAGE
	  DOUBLE PRECISION, SAVE, DIMENSION(:), POINTER :: Sndfw, Sodfw 
	  DOUBLE PRECISION, SAVE, DIMENSION(:), POINTER :: dGWQ, Q, leak 
	  DOUBLE PRECISION, SAVE, DIMENSION(:), POINTER :: Dfdh, KSTRM
	  DOUBLE PRECISION, SAVE, DIMENSION(:), POINTER :: KSB, BRHS,HCOEF
	  DOUBLE PRECISION, SAVE, POINTER :: Thickfactdfw
        DOUBLE PRECISION,SAVE,DIMENSION(:),  POINTER:: dqdhstr,WETTIME
        DOUBLE PRECISION,SAVE,DIMENSION(:),  POINTER:: SUMLEAK,SUMRCH
        DOUBLE PRECISION,SAVE,DIMENSION(:),  POINTER:: fderiv, Func
        DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER:: DELSTOR, WETPER
        DOUBLE PRECISION,SAVE,DIMENSION(:,:),POINTER:: STRMdfw,SEGDFW
	  DOUBLE PRECISION, SAVE, POINTER :: RMS1dfw, RMS2dfw, fheaddfw 
	  DOUBLE PRECISION, SAVE, POINTER :: Fheadsavedfw,ffluxdfw,RMSdfw
        INTEGER,SAVE,DIMENSION(:,:),POINTER:: ICELLDFW
	  INTEGER,SAVE,DIMENSION(:,:,:),POINTER:: IGWINDEX
      TYPE GWFDFWTYPE
        INTEGER, POINTER:: NSS, ISTCB1, ISTCB2
        INTEGER, POINTER:: MAXPTS, IRTFLG
        INTEGER, POINTER:: NSTRMdfw, NUMACTIVEDFW
        INTEGER, POINTER:: ITMP, IRDFLG, IPTFLG
        REAL, POINTER:: CONST
        INTEGER,       DIMENSION(:),  POINTER:: IOTSG, NSEGCK,WETTED
        INTEGER,       DIMENSION(:,:),POINTER:: ISEGDFW, IDIVA, ISTRMdfw
        REAL,          DIMENSION(:),  POINTER:: SGOTFL
        DOUBLE PRECISION,     DIMENSION(:,:),POINTER:: QSTAGE
	  DOUBLE PRECISION,     DIMENSION(:), POINTER :: Sndfw, Sodfw
	  DOUBLE PRECISION,     DIMENSION(:), POINTER :: dGWQ, Q, leak
	  DOUBLE PRECISION,     DIMENSION(:), POINTER :: Dfdh, KSTRM
        DOUBLE PRECISION,     DIMENSION(:), POINTER :: KSB, BRHS, HCOEF
	  DOUBLE PRECISION,     POINTER :: Thickfactdfw
        DOUBLE PRECISION,     DIMENSION(:),  POINTER:: dqdhstr, WETTIME
        DOUBLE PRECISION,     DIMENSION(:),  POINTER:: SUMLEAK, SUMRCH
        DOUBLE PRECISION,     DIMENSION(:),  POINTER:: fderiv, Func
        DOUBLE PRECISION,     DIMENSION(:,:),POINTER:: DELSTOR, WETPER
        DOUBLE PRECISION,     DIMENSION(:,:),POINTER:: STRMdfw, SEGDFW
        DOUBLE PRECISION, POINTER :: RMS1dfw, RMS2dfw, RMSdfw
	  DOUBLE PRECISION, POINTER :: Fheadsavedfw, ffluxdfw, fheaddfw
	  INTEGER,DIMENSION(:,:),POINTER:: ICELLDFW
	  INTEGER,DIMENSION(:,:,:),POINTER:: IGWINDEX
      END TYPE
      TYPE(GWFDFWTYPE), SAVE:: GWFDFWDAT(10)
      END MODULE GWFDFWMODULE
C
C-------SUBROUTINE GWF2DFW1AR
      SUBROUTINE GWF2DFW1AR(In, Iunitbcf, Iunitlpf, Iunithuf,  
     +                      Iunitnwt, Iouts, Igrid)   
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR STREAMS
C     INITIALIZE VARIABLES FOR DFW PACKAGE
C     READ STREAM DATA THAT IS CONSTANT FOR ENTIRE SIMULATION:
C     REACH DATA 
C     VERSION  1.0.00: 2009
C     ******************************************************************
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFDFWMODULE
      USE GLOBAL,       ONLY: IOUT, IBOUND, BOTM, STRT, DELR, DELC,
     +                        ITRSS, NLAY, NROW, NCOL
      USE GWFLPFMODULE, ONLY: SC2LPF=>SC2
      USE GWFBCFMODULE, ONLY: SC1, SC2, LAYCON
      USE GWFHUFMODULE, ONLY: SC2HUF
	USE GWFNWTMODULE, ONLY: IA, JA, A
      IMPLICIT NONE
      INTRINSIC ABS, DBLE
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER In, Iunitbcf, Iunitlpf, Iunithuf, Iouts, 
     +        Igrid, Iunitnwt
C     ------------------------------------------------------------------
C     LOCAL VARIABLES
C     ------------------------------------------------------------------
      CHARACTER*200 line
      INTEGER lloc, istart, istop, i, ii, nlst, lb, ichk, icalc, ij, jj
      INTEGER nseg, nreach, krch, irch, jrch, jseg, ireach, ksfropt
      INTEGER krck, irck, jrck, jsegck, ireachck, kkptflg, ib,nper
      INTEGER lstsum, lstbeg, numinst, idum(1), ip, iterp, mstrmar
      INTEGER nssar, nstrmar, nsegdim, itrib, irchtrib, idiv, irchdiv
      REAL r, seglen, sumlen, thsslpe, thislpe, uhcslpe, rchlen, dist
      REAL epsslpe, runoff, etsw, pptsw
C     ------------------------------------------------------------------
      iterp = 1
      idum(1) = 0
      ALLOCATE (NSS, TOTSPFLOW, HEPS, CLOSEZERO, QTOL)
      ALLOCATE (ISTCB1, ISTCB2, MAXPTS)
      ALLOCATE (NSTRMdfw,Numactivedfw)
      ALLOCATE (ITMP, IRDFLG, IPTFLG)
      ALLOCATE (CONST, IRTFLG)
      ALLOCATE (THICKFACTDFW)
      ALLOCATE (RMS1dfw, RMS2dfw, RMSdfw)
      RMS1dfw = 0.0D0
      RMSdfw = 0.0D0
	THICKFACTDFW = 0.0001
      HEPS = 1.0E-7
	Qtol = 1.0E-4
      CLOSEZERO = 1.0E-15

C1------IDENTIFY PACKAGE AND INITIALIZE NSTRMdfw.
      WRITE (IOUT, 9001) In
 9001 FORMAT (1X, /, ' DFW1 -- DIFFUSION WAVE PACKAGE, '
     +        ,'VERSION 1, 08/15/2010', /, 9X, 
     +         'INPUT READ FROM UNIT', I4)
C
C2------READ COMMENT RECORDS, NSTRMdfw, NSS, DLEAK, ISTCB1, ISTCB2.

      CALL URDCOM(In, IOUT, line)
      lloc = 1
      IRTFLG = 0
      CALL URWORD(line, lloc, istart, istop, 2, NSTRMdfw, r, IOUT, In)
      CALL URWORD(line, lloc, istart, istop, 2, NSS, r, IOUT, In)
      CALL URWORD(line, lloc, istart, istop, 3, i, CONST, IOUT, In)
      CALL URWORD(line, lloc, istart, istop, 2, ISTCB1, r, IOUT, In)
      CALL URWORD(line, lloc, istart, istop, 2, ISTCB2, r, IOUT, In)
      
      ALLOCATE (Icelldfw(NSTRMdfw,NSS),Func(NSTRMdfw))
	ALLOCATE (Igwindex(NCOL, NROW, NLAY))
	ALLOCATE (JA(NSS),IA(NSTRMdfw*7))
	ALLOCATE (Fderiv(NSTRMdfw*7),Kstrm(NSTRMdfw*7))
	ALLOCATE (Sndfw(NSTRMdfw),Sodfw(NSTRMdfw),Dfdh(NSTRMdfw))
	ALLOCATE (Ksb(NSTRMdfw), BRHS(NSTRMdfw), HCOEF(NSTRMdfw))
	ALLOCATE (Fheaddfw, Ffluxdfw, Fheadsavedfw) 
	ALLOCATE (dGWQ(NSTRMdfw), dqdhstr(NSTRMdfw), Q(NSTRMdfw))
	ALLOCATE (LEAK(NSTRMdfw))
	Q = 0.0d0
	LEAK = 0.0D0  
      Icelldfw = 0
	igwindex = 0
	Sndfw = 0.0D0
      Sodfw = 0.0D0
      Dfdh = 0.0D0
	Kstrm = 0.0D0
	BRHS = 0.0D0
	HCOEF = 0.0D0
	Ksb = 0.0
      Fheaddfw = 0.0D0
      Fheadsavedfw = 0.0D0
	Func = 0.0D0
	dGWQ = 0.0D0
	dqdhstr = 0.0D0
	Numactivedfw = NSTRMdfw
!
!4b-----LEFTOVER FROM SFR
!      CALL URWORD(line, lloc, istart, istop, 2, IRTFLG, r, IOUT, In)
      IF ( NSS.LT.0 ) NSS = 0
      nssar = 1
      IF (NSS.GT.0) nssar = NSS
      nstrmar = 1
      IF (NSTRMdfw.GT.0) nstrmar = NSTRMdfw
      nsegdim = NSS 
      IF (nsegdim.LT.1) nsegdim = 1
C
C5------CALCULATE SPACE NEEDED FOR TABULATED DISCHARGE VERSUS FLOW
C         AND WIDTH RELATIONS.
      MAXPTS = 5*20
C
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR STREAMS
C     ******************************************************************
C
      ALLOCATE (STRMdfw(34,nstrmar), ISTRMdfw(5,nstrmar))
	ALLOCATE (WETTIME(nstrmar),WETTED(nstrmar))
	WETTED = 0
	WETTIME = 0.0
      STRMdfw = 0.0
      ISTRMdfw = 0
      ALLOCATE (SEGDFW(26,nsegdim),ISEGDFW(4,nsegdim),IDIVA(2,nsegdim))
      SEGDFW = 0.0D0
      ISEGDFW = 0
      IDIVA = 0
      ALLOCATE (IOTSG(nsegdim))
      IOTSG = 0
C
C6------PRINT INFORMATION THAT WAS READ.
      WRITE (IOUT, 9003) NSTRMdfw, NSS, CONST
!      IF ( ISTCB1.GT.0 ) WRITE (IOUT, 9006) ISTCB1
!      IF ( ISTCB2.GT.0 ) WRITE (IOUT, 9007) ISTCB2
 9003 FORMAT (//, ' NUMBER OF STREAM REACHES IS', I5, //, 
     +        ' NUMBER OF STREAM SEGMENTS IS', I5, //, 
     +        ' MAXIMUM ERROR FOR STREAM LEAKAGE RATES IS', 
     +        1PE10.2,///)
! 9006 FORMAT (' FLOW TO AND FROM GROUND WATER FOR EACH STREAM REACH ',
!     +        'WILL BE SAVED ON UNIT', I3)
! 9007 FORMAT (' STREAM OUTPUT WILL BE WRITTEN TO FILE ON UNIT', I4)
C
C7------CHECK FOR ERRORS.
      IF ( NSTRMdfw.LE.0 .OR. NSS.LE.0 ) THEN
        WRITE (IOUT, 9008)
        In = 0
        NSS = 0
        NSTRMdfw = 0
        RETURN
      END IF
 9008 FORMAT (//, 'NO STREAM REACHES (NSTRMdfw) AND/OR SEGMENTS (NSS)-',
     +        //, 'DFW PACKAGE BEING TURNED OFF'///)

      ALLOCATE (QSTAGE(MAXPTS,nsegdim)) 
      QSTAGE = 0.0
      ALLOCATE (NSEGCK(nssar), SGOTFL(nssar))
      NSEGCK = 0
      SGOTFL = 0.0
C
C10-----READ AND PRINT DATA FOR EACH STREAM REACH. 
      WRITE (IOUT, 9014)
 9014 FORMAT (1X, //, 3X, 'STREAM NETWORK DESCRIPTION: ', //, 3X,
     +        'LAYER    ROW    COL   SEGMENT   REACH     LENGTH',
     +        '     STREAMBED     STREAMBED   STREAMBED     STREAMBED',
     +        '    THETA S     THETA I    THETA R'
     +        /, 26X, 'NUMBER   NUMBER    IN CELL    TOP ELEV.    ', 
     +        '   SLOPE     THICKNESS', '   HYDR. CONDUCT.', /, 3X,
     +        105('-'))
C
C11-----READ AND WRITE DATA FOR EACH REACH 
      nseg = 0
      nreach = 0    
      DO ii = 1, NSTRMdfw
          READ (In, *) krch, irch, jrch, jseg, ireach, STRMdfw(1, ii), 
     +                 STRMdfw(3, ii), STRMdfw(2, ii), STRMdfw(8, ii), 
     +                 STRMdfw(6, ii), STRMdfw(32, ii), STRMdfw(33,ii),
     +		              STRMdfw(34, ii)
          STRMdfw(4, ii) = STRMdfw(3, ii) - STRMdfw(8, ii)
          IF ( STRMdfw(2, ii).LE.0.0 ) THEN
            WRITE (IOUT, 9017) jseg, ireach
            CALL USTOP(' ')
          END IF
 9017   FORMAT (//, ' ***ERROR***  SLOPE IS SPECIFIED LESS THAN OR ', 
     +          'EQUAL TO ZERO FOR SEGMENT', I8, ' REACH', I8, /, 
     +          ' PROGRAM IS STOPPING')
		 WRITE (IOUT, 9020) krch, irch, jrch, jseg, ireach, STRMdfw(1, ii),
     +                 STRMdfw(3, ii), STRMdfw(2, ii), STRMdfw(8, ii),
     +                 STRMdfw(6, ii), STRMdfw(32, ii), STRMdfw(33,ii), 
     +                 STRMdfw(34,ii)
 9020   FORMAT (2X, I6, 2I7, I8, I9, 3X, 1PE11.4, 2X, 1PE11.4, 2X, 
     +          1PE11.4, 2X, 1PE11.4, 4(2X, 1PE11.4))
C
C13-----CHECK RANGE AND ORDER FOR SEGMENTS AND REACHES.
        IF ( jseg.LE.0 .OR. jseg.GT.NSS ) THEN
          WRITE (IOUT, 9023)
          CALL USTOP(' ')
        END IF
        IF ( jseg.NE.nseg ) THEN
          nseg = nseg + 1
          nreach = 0
          IF ( jseg.NE.nseg ) THEN
            WRITE (IOUT, 9024)
            CALL USTOP(' ')
          END IF
        END IF
        nreach = nreach + 1
        IF ( ireach.NE.nreach ) THEN
          WRITE (IOUT, 9025)
          CALL USTOP(' ')
        END IF
 9023   FORMAT (' SEGMENT MUST BE GREATER THAN 0 AND LESS THAN NSS')
 9024   FORMAT (' SEGMENTS MUST BE IN ORDER FROM 1 THROUGH NSS')
 9025   FORMAT (' EACH SEGMENT MUST START WITH REACH 1, AND', /, 
     +          ' REACHES MUST BE NUMBERED CONSECUTIVELY')
        ISTRMdfw(1, ii) = krch
        ISTRMdfw(2, ii) = irch
        ISTRMdfw(3, ii) = jrch
        ISTRMdfw(4, ii) = jseg
        ISTRMdfw(5, ii) = ireach
        icelldfw(ireach, jseg) = ii
		   igwindex(jrch,irch,krch) = ii
        SEGDFW(1,ISTRMdfw(4,ii))=SEGDFW(1,ISTRMdfw(4,ii))+STRMdfw(1, ii)
C       Number of reaches in segment added to ISEG
        ISEGDFW(4, jseg) = ireach
        STRMdfw(5, ii) = 1.0
        STRMdfw(7, ii) = 1.0
!              STRMdfw(13, ii) = etsw*rchlen
!              STRMdfw(14, ii) = pptsw*rchlen
        STRMdfw(15, ii) = STRMdfw(3, ii)  !initial value of strmhd
!		         END IF
	  STRMdfw(31, ii) = STRMdfw(15,ii)    !set the old head too
      END DO
!
C18-----COMPUTE STREAM REACH VARIABLES.
      irch = 1
      ksfropt = 0
      DO nseg = 1, NSS
        seglen = SEGDFW(1, nseg)
        runoff = SEGDFW(3, nseg)
        etsw = SEGDFW(4, nseg)
        pptsw = SEGDFW(5, nseg)
        sumlen = 0.0
	End do
C
C14-----READ SEGMENT INFORMATION FOR FIRST STRESS PERIOD.
      READ (In, *) ITMP, IRDFLG, IPTFLG
      nlst = NSS
      lb = 1
      ichk = 1
      CALL SGWF2DFW1RDSEG(nlst, lb, In, NSEGCK, NSS, ichk, 1)
C
C18-----CHECK IF STREAM REACH IS IN ACTIVE CELL.
      kkptflg = 0
      DO ichk = 1, NSTRMdfw
        krck = ISTRMdfw(1, ichk)
        irck = ISTRMdfw(2, ichk)
        jrck = ISTRMdfw(3, ichk)
        jsegck = ISTRMdfw(4, ichk)
        ireachck = ISTRMdfw(5, ichk)
        IF ( IBOUND(jrck, irck, krck).EQ.0 ) THEN
          kkptflg = kkptflg + 1
          IF ( kkptflg.EQ.1 ) WRITE (IOUT, 9029) jsegck, ireachck, 
     +                              IBOUND(jrck, irck, krck), krck,
     +                              irck, jrck
        ELSE IF ( IBOUND(jrck, irck, krck).LT.0 ) THEN
          WRITE (IOUT, 9030) jsegck, ireachck, IBOUND(jrck, irck, krck),
     +                       krck, irck, jrck
        END IF
      END DO
      IF ( kkptflg.EQ.1 ) THEN
        WRITE (IOUT, 9031)
      ELSE IF ( kkptflg.GT.1 ) THEN
        WRITE (IOUT, 9032) kkptflg
      END IF
C-----MAS 14/1/16 READ IN TIMINGS--------------------------------------- 
      OPEN (unit=4,FILE = "timings.in")
        READ(4,*)nper    
        ALLOCATE(TSEG(1:nper))
        ALLOCATE(TRCH(1:nper))
        DO i = 1,nper
          READ(4,*) tseg(i),trch(i)
        END DO
      CLOSE(4)
!     -----------------------------------------------------------------
C
 9029 FORMAT (/, ' *** WARNING *** FIRST OCCURRENCE WHERE A ', 
     +        'STREAM REACH IS ASSIGNED TO AN INACTIVE CELL IS SEGMENT',
     +        I5, ' REACH NO.', I5, /, '  IBOUND ARRAY VALUE IS', I5,
     +        ' AT LAYER', I5, '; ROW', I5, '; COLUMN', I5, '.')
 9030 FORMAT (/, ' *** WARNING *** STREAM SEGMENT', I5, ' REACH NO.',
     +        I5, ' IS CONNECTED TO A CONSTANT HEAD CELL.'/,
     +        ' IBOUND ARRAY VALUE IS', I5, ' AT ', 'LAYER', I5,
     +        '; ROW', I5, '; COLUMN', I5, '.', /,
     +        ' NO STREAM LEAKAGE WILL BE ALLOWED-- SUGGEST ', 
     +        'REMOVING STREAM REACH FROM CELL OR CHANGE CELL ', 
     +        'TO VARIABLE HEAD.', /)
 9031 FORMAT (/, ' *** WARNING *** ONLY 1 STREAM REACH WAS ', 
     +        'ASSIGNED TO A CELL WHERE THE IBOUND ARRAY WAS ZERO.', /,
     +        ' PROGRAM SEARCHES FOR UPPERMOST ACTIVE CELL IN VERTICAL',
     +        ' COLUMN,IF ALL CELLS ARE INACTIVE, STREAM LEAKAGE WILL',
     +        ' NOT BE ALLOWED. ', /)
 9032 FORMAT (/, ' *** WARNING *** A TOTAL OF', I6, 'STREAM REACHES ', 
     +        'WERE ASSIGNED TO CELLS WHERE THE IBOUND ARRAY WAS ZERO.',
     +        /, ' PROGRAM SEARCHES FOR UPPERMOST ACTIVE CELL IN',
     +        ' VERTICAL COLUMN FOR ALL OCCURRENCES.', /,
     +        ' IF ALL CELLS IN A VERTICAL COLUMN ARE INACTIVE,',
     +        ' STREAM LEAKAGE WILL NOT BE ALLOWED FOR ASSOCIATED',
     +        ' STREAM REACH. ', /)
C     ------------------------------------------------------------------
C23-----SAVE POINTERS FOR GRID AND RETURN.
      CALL SGWF2DFW1PSV(Igrid)
      RETURN
      END SUBROUTINE GWF2DFW1AR
C
C-------SUBROUTINE GWF2DFW1RP
      SUBROUTINE GWF2DFW1RP(In, Iunitlak, Kkper, Kkstp, Iouts, Iunitnwt,
     +                      Igrid)
C     ******************************************************************
C     READ STREAM DATA FOR STRESS PERIOD
C     VERSION  1.0.00: 2009
C     ******************************************************************
      USE GWFDFWMODULE  
      USE GLOBAL,       ONLY: IOUT, ISSFLG, IBOUND, BOTM, HNEW,NLAY
		 USE GWFNWTMODULE, ONLY: Numactive, IA,JA, Diag, Icell, Hiter
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER Kkper, Kkstp, In, Iunitlak, Iouts, Igrid, Iunitnwt
      INTEGER Ischk, Nlst
      DIMENSION Ischk(NSS)
C     ------------------------------------------------------------------
C     LOCAL VARIABLES
C     ------------------------------------------------------------------
      DOUBLE PRECISION h, sbot
      REAL avdpth, avthk, bottom, dist, dpslpe, dpth1, dpth2,
     +     dpthlw, dndiff, eldn, elslpe, wdth1, wdth2, wdthlw,
     +     hcslpe, rchlen, rough, roughbnk, roughch, 
     +     seglen, strlen, sumlen, thkslpe, top, wdslpe, arealw, 
     +     width, updiff, zero, hydr1, hydr2, hydrlw, area1,
     +     area2, mann1, mann2, mannlw
      INTEGER i, ic, icalc, ichk, icp, iflginit, ii, ik, il, ilay, ip,
     +        ipt, ir, irch, irp, isoptflg, iss, istep, istsg, iwvcnt,
     +        jj, jk, k5, k6, k7, kk, ksfropt, kss, ktot, l, lstbeg,
     +        nseg, nstrpts, idum, iqseg, iupseg, lstend, n, noutseg,
     +        nrea, ij, ireach,jseg, itrib, irchtrib, idiv,irchdiv, jl,
     +        jr, jc, iupsg, iuprch, starthere
C     ------------------------------------------------------------------
C  
C-------SET POINTERS FOR CURRENT GRID.
      CALL SGWF2DFW1PNT(Igrid)
C
C1------READ ITMP FLAG TO REUSE NON-PARAMETER DATA, 2 PRINTING FLAGS,
C         AND NUMBER OF PARAMETERS BEING USED IN CURRENT STRESS PERIOD. 
      iss = ISSFLG(Kkper)
      zero = 1.0E-7
	icalc = 4
      IF ( Kkper.GT.1 ) THEN
          READ (In, *) ITMP, IRDFLG, IPTFLG
      END IF
C
C6------READ STREAM SEGMENT DATA.
      IF ( ITMP.GE.0 ) THEN
        lstbeg = 1
        ichk = 1
        IF ( Kkper.GT.1 ) CALL SGWF2DFW1RDSEG(ITMP, lstbeg, In, 
     +                                        NSEGCK, NSS,
     +                                        ichk, Kkper)
      ELSE 
        WRITE (IOUT, 9003)
        RETURN
      END IF

 9001 FORMAT (/, ' CANNOT SPECIFY MORE THAN NSS STREAM SEGMENTS')
 9002 FORMAT (//, '  ***  STREAM SEGMENTS MUST BE DEFINED FOR ', 
     +       'FIRST STRESS PERIOD; CODE STOPPING ***') 
 9003 FORMAT (/, ' REUSING STREAM SEGMENT DATA FROM LAST STRESS PERIOD')
C
C8------CHECK FOR ERRORS IN SEGMENT DATA.
      IF ( ITMP.GT.0 .OR. kkper.EQ.1) THEN
        DO nseg = 1, NSS

C9------READ DATA ACCORDING TO VARIABLE ISFROPT.
!MAS isfropt deleted, always = 1
              isoptflg = 1
          IF ( IDIVA(2, nseg).GT.0 ) THEN
            WRITE (IOUT, 9008) nseg
            IDIVA(2, nseg) = 0
          ELSE IF ( IDIVA(2, nseg).LT.-3 ) THEN
            WRITE (IOUT, 9009) nseg
            IDIVA(2, nseg) = 0
          ELSE IF ( IDIVA(2, nseg).EQ.-2 ) THEN
            IF ( SEGDFW(2,nseg).LT.0.0.OR.SEGDFW(2,nseg).GT.1.0 ) THEN
              WRITE (IOUT, 9010) nseg
              SEGDFW(2, nseg) = 0.0
            END IF
          END IF
        END DO
 9008   FORMAT (/, 5X, '*** WARNING *** IPRIOR > 0 FOR NSEG =', I7, /,
     +          10X, 'THIS OPTION NOT YET AVAILABLE; CODE WILL ', 
     +          'ASSUME IPRIOR = 0', /)
 9009   FORMAT (/, 5X, '*** WARNING *** IPRIOR < -3 FOR NSEG =', I7, /,
     +          10X, 'THIS VALUE IS OUT OF RANGE; CODE WILL ', 
     +          'ASSUME IPRIOR = 0', /)
 9010   FORMAT (/, 5X, '*** WARNING *** IPRIOR = -2 FOR NSEG =', I7,
     +          ' & FLOW VALUE IS OUT OF RANGE (.0 - 1.);', /, 10X,
     +          'ASSUME NO DIVERSION OF FLOW', /)
C
C10-----PLACE STREAM SEGMENT IDENTITY NUMBERS IN ISEG ARRAY.
C         5 ASSIGNED TO SEGMENTS NOT RECEIVING TRIBUTARY FLOW.
C         6 ASSINGED TO SEGMENTS THAT DIVERT FLOW.
C         7 ASSIGNED TO SEGMENTS RECEIVING TRIBUTARY FLOW.
        k5 = 0
        k6 = 0
        k7 = 0
        DO nseg = 1, NSS
C
C11-----IDENTIFY SEGMENTS THAT DIVERT FLOW.
          IF ( IDIVA(1, nseg).NE.0 ) THEN
            ISEGDFW(3, nseg) = 6
            k6 = k6 + 1
C
C12-----IDENTIFY SEGMENTS THAT DO NOT DIVERT FLOW.
          ELSE
            jj = 0
C
C13-----IDENTIFY SEGMENTS THAT RECEIVE TRIBUTARY FLOW.
            DO ii = 1, NSS
              IF ( IOTSG(ii).EQ.nseg ) jj = 1
            END DO
C
C14-----IDENTIFY SEGMENTS THAT DO NOT RECEIVE TRIBUTARY FLOW.
            IF ( jj.EQ.0 ) THEN
              ISEGDFW(3, nseg) = 5
              k5 = k5 + 1
            ELSE
              ISEGDFW(3, nseg) = 7
              k7 = k7 + 1
              IF ( jj.NE.1 ) WRITE (IOUT, 9011) nseg, jj
            END IF
          END IF
        END DO
C
C15-----TALLY DIFFERENT STREAM SEGMENT TYPES - this is for accounting only.
        ktot = k5 + k6 + k7
        WRITE (IOUT, 9012) k5, k6, k7
C
C16-----PRINT WARNING IF TALLIED SEGMENTS LESS THAN NSS.
        IF ( ktot.NE.NSS ) THEN
          WRITE (IOUT, 9013) ktot, NSS
          CALL USTOP(' ')
        END IF
 9011   FORMAT (//, 5X, '*** WARNING *** ERROR WHILE ', 
     +          'CLASSIFYING SEGMENTS:   NSEG =', I6, 4X, 'JJ =', I6,//)
 9012   FORMAT (///1X, 'CLASSIFICATION & COUNT OF STREAM SEGMENTS ', 
     +          'BASED ON SOURCE OF INFLOW:', //, 16X, 
     +          'HEADWATER     DIVERSION     RECEIVES TRIBUTARY FLOW', /
     +          16X, '---------     ---------    ', 
     +          ' -----------------------', /, 16X, I6, I15, I16, /)
 9013   FORMAT (/, 5X, '*** WARNING ***  INTERNAL ERROR SUMMING ', 
     +          'TYPES OF STREAM SEGMENTS:  NSEG =', I6, 5X, 'JJ =',
     +          I6//)
C
C17-----PRINT INPUT DATA IF IRDFLG IS ZERO.
        IF ( IRDFLG.LE.0 ) CALL SGWF2DFW1PRSEG(NSS, 1, Kkper,Iouts)

      END IF
C
! Fill CRS pointers for MODFLOW storage scheme by relating Row, Col, and Lay to
! Jacobian order
! Modified from NWT fillindex
!ij is the number of active cells (row in sol. vector)
!jj is the number of non-zero elements in the Jacobian
!IA() is the pointer to a new row in Jacobian (CRS)
!JA() is the order in the unknown vector (CRS)
      IF ( Kkper.EQ.1 ) THEN
        jj = 1
        DO ij = 1, NSTRMdfw
		     jl = ISTRMdfw(1, ij)
          jr = ISTRMdfw(2, ij)
          jc = ISTRMdfw(3, ij)
          ireach = ISTRMdfw(5, ij)
          jseg = ISTRMdfw(4, ij)
! this is not part of fillindex; just using do-loop to set inital hiter
		Hiter(ij) = STRMdfw(3, ij)
          IA(ij) = jj
		JA(jj) = Icelldfw(ireach, jseg)
          jj = jj + 1
 ! connected to upstream reach of same segment
          IF ( ireach.GT.1 ) THEN
            JA(jj) = Icelldfw(ireach-1,jseg)
            jj = jj + 1
! diversion connected to upstream segment
		ELSE IF (IDIVA(1, jseg).GT.0) THEN
		  iupsg = IDIVA(1,jseg)
		  iuprch = ISEGDFW(4, iupsg)
            JA(jj) = Icelldfw(iuprch, iupsg)
            jj = jj + 1
! connected to tributary segments.
          ELSE
		  DO itrib = 1, NSS
              IF ( jseg.EQ.IOTSG(itrib) ) THEN
                irchtrib = ISEGDFW(4,itrib)
                JA(jj) = Icelldfw(irchtrib, itrib)
                jj = jj + 1
		    END IF
		  END DO   
          END IF
 ! connected to outflowing diversion
          IF ( ireach.EQ.ISEGDFW(4, jseg))THEN 
            DO idiv = 1, NSS		 		 		 		 
		     IF ( jseg.EQ.IDIVA(1, idiv) ) THEN 
		       irchdiv = 1                      
		       JA(jj) = Icelldfw(irchdiv, idiv)
		       jj = jj + 1
		     END IF
		  END DO		 
		END IF   
 ! Test for downstream reaches, is it the last reach in a segment         
          IF (ireach.LT.ISEGDFW(4, jseg) )  THEN
		  JA(jj) = Icelldfw(ireach + 1, jseg)
		  jj = jj + 1
		ELSEIF ( IOTSG(jseg).GT.0 ) THEN
		  JA(jj) = Icelldfw(1, IOTSG(jseg))
		  jj = jj + 1
		ENDIF
 ! connected to groundwater; this is just a place holder, JA(jj) from NWT
          IF ( IBOUND(jc, jr, jl).GT.0 ) THEN
		  JA(jj) = Icelldfw(ireach, jseg)    
		  jj = jj + 1
		END IF
        ENDDO
        IA(NSTRMdfw+1) = jj 
      END IF
      RETURN
      END SUBROUTINE GWF2DFW1RP
C
C-------SUBROUTINE GWF2DFW1FM
      SUBROUTINE GWF2DFW1FM(Kkiter,Kkper,Kkstp,Iunitlak,Iunitnwt,Igrid)
C     *****************************************************************
C     VERSION  1.0.00:  2010
C     Builds Jacobian and then calls NWT to solve it
C     *****************************************************************
      USE GWFDFWMODULE
      USE GLOBAL,       ONLY: NLAY, IOUT, ISSFLG, IBOUND, HNEW
      USE GWFBASMODULE, ONLY: DELT, TOTIM, HDRY
!      USE GWFLAKMODULE, ONLY: Nlakesar, Theta, Stgold, Stgnew, Vol
      USE GWFNWTMODULE, ONLY: IA,JA, Hchange, A, BB
      IMPLICIT NONE
      INTRINSIC IABS, ABS, DABS, MIN, DSQRT, FLOAT, SQRT, SNGL
C     -----------------------------------------------------------------
C     ARGUMENTS
C     -----------------------------------------------------------------
      INTEGER Kkiter, Kkper, Iunitlak, Igrid, Kkstp, Iunitnwt
C     -----------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
      DOUBLE PRECISION, EXTERNAL ::CONFUN,dfwB,SBCON 
C     LOCAL VARIABLES
C     -----------------------------------------------------------------
      INTEGER irchhld, iseghld, ireach, jseg, ireachcon, isegcon
      INTEGER ij, itstr, l, lfold,lsub, ll, il, ir, ic, iss, kflag,
     +        nstrpts, jj, istsg, nreach, ijnow, ijup, ijcon, ik
	DOUBLE PRECISION dvrsn, slope, Hcon, Hx, Hup,Hold, Hc, Hr
      DOUBLE PRECISION Hnow, Hsb, ferr
C     -----------------------------------------------------------------
C
C-------SET POINTERS FOR CURRENT GRID.
      CALL SGWF2DFW1PNT(Igrid)
      ISS = ISSFLG(Kkper)  
	IF (kkper == 62) THEN
	 kflag = 0
	END IF  
C
C1------RETURN IF NO STREAMS (NSTRMdfw<=0).
      IF ( NSTRMdfw.LE.0 ) RETURN
      irchhld = 1
      iseghld = 1
      Ffluxdfw = 0.0D0
      Fheadsavedfw = Fheaddfw
      Fheaddfw = 0.0D0
      kflag = 0
      CALL Sndfw_update() 
	DO jj = 1, NSTRMdfw
	  IF(iss ==0.AND.kkiter==1.AND.kkper.GT.1) THEN
		 IF (WETTED(jj)==0) THEN
		   WETTIME(jj)= delt
		 ELSE
		   WETTIME(jj) = WETTIME(jj) + delt
		 END IF
	  END IF
	  Hnow = STRMdfw(15, jj) 
	  Hold = STRMdfw(31, jj) 
	  DO ik = IA(jj),IA(jj+1)-1
		A(ik) = 0.0D0
        End Do
	  BB(jj) = 0.0D0
	  IF (ISS == 0) Then    
          CALL STORFM(Hr, Hc, Hnow, Hold, jj)
!		       Hcoef(jj)= Hc
!		       Brhs(jj) =  Hr
	  END IF 
        ireach = ISTRMdfw(5, jj)
        jseg = ISTRMdfw(4, jj)         
	  il = ISTRMdfw(1, jj)
        ir = ISTRMdfw(2, jj)
        ic = ISTRMdfw(3, jj)		 
        STRMdfw(19,jj) = HNEW(ic,ir,il)		   
	  Hsb = STRMdfw(19, jj)
	  If (Hnow .GT. Hsb) THEN
		Hup = Hnow
	  ELSE
		Hup = Hsb
	  END IF
	  Ksb(jj) = SBCON(jj, Hnow, Hsb, Hup,iss)      
! Check all connectivities in the reach, get conductances
	  DO ij = IA(jj)+1,IA(jj+1)-2 
		 ireachcon = ISTRMdfw(5, JA(ij))
		 isegcon = ISTRMdfw(4, JA(ij))
		 Hcon = STRMdfw(15, JA(ij)) 
		 IF ( Hcon-Hnow.GT.CLOSEZERO)THEN
             ijup = JA(ij)
             Kstrm(ij) = CONFUN(ijup, jj, JA(ij), Hcon, kflag)
		 ELSE		     
             ijup = jj
             Kstrm(ij) = CONFUN(ijup, jj, JA(ij), Hnow, kflag)
           ENDIF
        END DO
	  Func(jj) = dfwB(jj, Hnow, iss) - (Hc*Hnow-Hr)
! Find biggest headchange
	  IF ( ABS(Hchange(jj)).GT.ABS(fheaddfw) ) THEN
           Fheaddfw = Hchange(jj)
           irchhld = ISTRMdfw(5, jj)
           iseghld = ISTRMdfw(4, jj) 
        END IF 
      END DO
      RMS2dfw = RMS1dfw 
      CALL JACOBYAN(ISS)
      OPEN(991,file='DFW_solve')
      WRITE(991,111)
  111 FORMAT (1X,'Reach   Seg  Stper  TmStp  Itertn',
     +          '       RMS-New               RMS-Old           ',
     +          'Solver-Max-Delh' )
      WRITE (991, 9001) irchhld,iseghld,kkper,kkstp,kkiter,RMS1dfw,
     +                          RMS2dfw,FHEADSAVEDFW  
 9001 FORMAT (1x,I4,3X,I2,3(3X,I4), E20.10,4(2X,E20.10))
      END SUBROUTINE GWF2DFW1FM
!----------------------------------------------------------------------
      SUBROUTINE GWF2DFW1BD(Kkiter,Kkstp,Kkper,Iunitgage,Igrid)
C     *****************************************************************
C     CALCULATE VOLUMETRIC GROUND-WATER BUDGET FOR STREAMS AND SUM
C     STREAMFLOWS IN MODELED AREA
C     VERSION  7.1.01: February 15, 2009
C     *****************************************************************
      USE GWFDFWMODULE
      USE GLOBAL,       ONLY: NCOL, NROW, NLAY, IOUT, ISSFLG, IBOUND,
     +                        HNEW, BUFF, NSTP
      USE GWFBASMODULE, ONLY: MSUM, ICBCFL, IBUDFL, DELT, PERTIM, TOTIM,
     +                        VBVL, VBNM, HDRY
	USE GWFNWTMODULE, ONLY: BB, JA, IA
      IMPLICIT NONE
      INTRINSIC FLOAT, ABS, IABS, DSQRT, DLOG10, SQRT, SNGL
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
C     FUNCTIONS
C     ------------------------------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: dfwB, SBCON, LOOKUP, CONFUN
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER Kkstp, Kkper, Iunitgage, maxiter
      INTEGER Igrid, ICNVG, Kkiter
C     ------------------------------------------------------------------
C     LOCAL VARIABLES
C     ------------------------------------------------------------------
      REAL zero, rtime
      INTEGER i, ibd, iblst, ibdlbl, ibdst, ibstlb, ic, ij, jj, 
     +        il, ilay, iout1, iout2, iprior, iprvsg, ir, istsg, itrib,
     +        kss, l, lk, ll, nreach, numdelt, iuprch, iupsg,pern, 
     +        nstrpts, iss, lsub, irt, itstr, IDIM, idiv, irchdiv,
     +        ireachcon, isegcon, ijup, irchtrib, ijcon, kflag
      DOUBLE PRECISION flowin,flobot, flow, flowot, fkstrm, qdiv, qin,
     +                 width, runof, runoff, precip, etstr, qout, AREA,
     +                 slope, depth, totflwt, totdelstor, hold, 
     +                 thetas, epsilon, thr, hnow,hup, hcon, Hr, Hc,
     +	                 Areacum, Leakcum, Perimcum, P
      DOUBLE PRECISION Hx, Hsb, Hnow1, STOR1, RATIN,RATOUT,nearzero
C     ------------------------------------------------------------------
C     LOCAL STATIC VARIABLES
C     ------------------------------------------------------------------
      CHARACTER*16 text, strtxt, txtlst
      DATA text/'  STREAM LEAKAGE'/
      DATA strtxt/'STREAMFLOW OUT  '/
      DATA txtlst/'STREAM LISTING  '/
C     -----------------------------------------------------------------
C-------SET POINTERS FOR THE CURRENT GRID.
      CALL SGWF2DFW1PNT(Igrid)
      iss = ISSFLG(Kkper)
	Areacum = 0.0
	Leakcum = 0.0
      ibd = 0
      ibdst = 0
      iblst = 0
      stor1 = 0.0D0
      kflag = 0
	NearZERO = 1.0E-09
	ZERO = 0.0
	RATOUT = 0.0
      RATIN = 0.0
C2------WRITE HEADER WHEN CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST.
      IF ( ibd.EQ.2 ) CALL UBDSV2(Kkstp, Kkper, text, iout1, NCOL, NROW,
     +                            NLAY, NSTRMdfw, IOUT, DELT, PERTIM, 
     +                            TOTIM, IBOUND)
      IF ( ibdst.EQ.2 ) CALL UBDSV2 (Kkstp, Kkper, strtxt, iout2, NCOL, 
     +                               NROW, NLAY, NSTRMdfw, IOUT, DELT, 
     +                               PERTIM, TOTIM, IBOUND)
C
C3------CLEAR BUFFERS.
!      DO il = 1, NLAY
!        DO ir = 1, NROW
!           DO ic = 1, NCOL
!             BUFF(ic, ir, il) = zero
!           END DO
!        END DO
!      END DO
c
      OPEN(912,file='DFW_out')
      IF (kkstp == NSTP(KKPER)) THEN
	  WRITE (912, 9009) txtlst, Kkper,kkstp
      END IF
C5b------DETERMINE LAYER, ROW, COLUMN, SEG & RCH.
      DO l = 1, NSTRMdfw
	  ll = l - 1
        il = ISTRMdfw(1, l)
        ir = ISTRMdfw(2, l)
        ic = ISTRMdfw(3, l)
        istsg = ISTRMdfw(4, l)
        nreach = ISTRMdfw(5, l)
  	  STRMdfw(19,l) = HNEW(ic,ir,il)		   
	  Hsb = STRMdfw(19, l)
	  Hnow = STRMdfw(15, l)
	  Hold = STRMdfw(31, l)
	  fKstrm = 0.0
	  Qdiv = 0.0
	  CALL Sndfw_update() 
!Recalculate leakage  =volumetric Q of hchan-hgw
	  If (Hnow .GT. Hsb) THEN
		Hup = Hnow
	  ELSE
		Hup = Hsb
	  END IF
	  Ksb(l) = SBCON(l, Hnow, Hsb, Hup, iss)   
!Calculate Qin and Qout for each reach.
!Qin for first reach of segment is specified
	  IF (ISTRMdfw(5,l)==1)THEN
		 Q(l) = SEGDFW(2,ISTRMdfw(4, l)) 
!add in trib/US segmt inflow
		 DO itrib = 1, NSS
              IF ( istsg.EQ.IOTSG(itrib) ) THEN  !it has a trib
                irchtrib = ISEGDFW(4,itrib)      !iseg(4 = tot# reaches
                ijcon = Icelldfw(irchtrib, itrib) 
                Hcon = STRMdfw(15, ijcon) 
		      IF ( Hcon-Hnow.GT.CLOSEZERO)THEN
                  ijup = ijcon
                  fKstrm = CONFUN(ijup, l, ijcon, Hcon, kflag)
		      ELSE		     
                  ijup = l
                  fKstrm = CONFUN(ijup, l, ijcon, Hnow, kflag)
                ENDIF
		      Q(l) = Q(l) + fKstrm*sqrt(abs(Hcon-Hnow))
		    END IF
		 END DO
!add in inflow from inflowing diversions
		 IF (IDIVA(1, istsg).GT.0) THEN
		   iupsg = IDIVA(1,istsg)  
		   iuprch = ISEGDFW(4, iupsg)
             ijcon = Icelldfw(iuprch, iupsg)
	       Hcon = STRMdfw(15, ijcon)
		   IF ( Hcon-Hnow.GT.CLOSEZERO)THEN
               ijup = ijcon
               fKstrm = CONFUN(ijup, l, ijcon, Hcon, kflag)
		   ELSE		     
               ijup = l
               fKstrm = CONFUN(ijup, l, ijcon, Hnow, kflag)
             ENDIF
		   Q(l) = Q(l) + fKstrm*sqrt(abs(Hcon-Hnow))
		 END IF
!outflow for first reach
		 Hcon = STRMdfw(15, JA(IA(l+1)-2)) 
		 IF ( Hcon-Hnow.GT.CLOSEZERO)THEN
             ijup = JA(IA(l+1)-2)
             fKstrm = CONFUN(ijup,l,JA(IA(l+1)-2),Hcon,kflag)
		 ELSE		     
             ijup = l
             fKstrm = CONFUN(ijup,l,JA(IA(l+1)-2),Hnow,kflag)
           ENDIF
		 Qout=fKstrm*sqrt(abs(Hcon-Hnow))
!else not first reach in segment
        ELSE
		 DO ij = IA(l)+1,IA(l+1)-2 
		    Hcon = STRMdfw(15, JA(ij)) 
		    IF ( Hcon-Hnow.GT.CLOSEZERO)THEN
                ijup = JA(ij)
                Kstrm(ij) = CONFUN(ijup, l, JA(ij), Hcon, kflag)
		    ELSE		     
                ijup = l
                Kstrm(ij) = CONFUN(ijup, l, JA(ij), Hnow, kflag)
              ENDIF
	!upstream
              IF (ij == IA(l)+1) THEN
			    Q(l)=Kstrm(ij)*sqrt(abs(Hcon-Hnow))
	!downstream - should be the same as Q(l) - Ksb 
		      ELSEIF (nreach.LT.ISEGDFW(4,istsg).AND.ij==IA(l+1)-2) THEN 
      		    Qout=Kstrm(ij)*sqrt(abs(Hcon-Hnow))
	!only other option is a diversion    !MAS030916 next 3 lines were screwing up Qout for last rch of ea seg
	       ELSEIF (nreach==ISEGDFW(4,istsg).AND.IDIVA(1, istsg).GT.0) THEN  !MAS030916 added in IF statements- does this properly capture outfl divsns?
	          Qdiv = Kstrm(ij)*sqrt(abs(Hcon-Hnow))
	        END IF
		 END DO
        END IF
! Set up the rest of the values for the .out file 
C8------STORE OUTFLOW FROM PREVIOUS SEGMENT IN SGOTFLW LIST FOR GAGES - gage doesn't work yet
        IF ( istsg.GT.1 ) THEN
           iprvsg = ISTRMdfw(4, ll)
           SGOTFL(iprvsg) = STRMdfw(9, ll)
	  END IF
C
C26-----DETERMINE STREAM WIDTH AND DEPTH.
        depth = Hnow - STRMdfw(3, l)
	  nstrpts = ISEGDFW(2, istsg)      
	  IDIM = nstrpts+1     !B goes from nstrpts+1 to 2*nstrpts
	  Hnow1 = Hnow -  STRMdfw(3,l) 
	  width = LOOKUP (Hnow1, istsg, IDIM, nstrpts)
	  IDIM = 3*nstrpts +1
	  area =  LOOKUP (Hnow1, istsg, IDIM, nstrpts)
	  STRMdfw(5, l) = width
	  IDIM = 4*nstrpts+1        
        P =  LOOKUP (Hnow1, istsg, IDIM, nstrpts)
        runof = 0.0D0       !STRM(12, l)
        runoff = 0.0D0      !STRM(24, l)
        precip = STRMdfw(14, l)*width
        etstr = STRMdfw(13, l)*width
C2------IF THE STRESS PERIOD IS TRANSIENT, ADD STORAGE TO HCOF AND RHS
        IF(ISS.EQ.0) THEN
          CALL STORBD(Hr, Hc, Hnow, Hold, l)
          Hcoef(l)= Hc
		Brhs(l) =  Hr
          Stor1 = Hcoef(l) - Brhs(l)
        END IF
! out for last reach
	  IF (nreach == ISEGDFW(4,istsg)) THEN
         Qout=Q(l) + Ksb(l) - stor1 -Qdiv !MAS 9/8/10 removed abs MAS030916 
         !write(912,*) Qout, Ksb(l), Qdiv, stor1
        end If
        STRMdfw(9,l) = Qout
! SAVE PREVIOUS ITERATION HEADS and Saturations        
        STRMdfw(31, l) = STRMdfw(15, l)
        Sodfw(l) = Sndfw(l)
! Set inflow and outflow less than some number to zero in printout
	  IF(Qout.LT.QTOL) Qout = 0.0
	  IF(Q(l).LT.QTOL) Q(l) = 0.0
 !     IF(STOR1.LT.QTOL) STOR1 = 0.0
	  IF (Q(l).GT.Nearzero) WETTED(l)=1
!MAS 13/1/15 MOVED NEXT 7 LINES ABOVE THE PRINT STATEMENT- DOES THIS FIX NEG Q?		   
        IF(Ksb(l).LE.ZERO) THEN
C5G-----FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
          RATIN = RATIN - Ksb(l)
        ELSE
C5H-----FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT. 
          RATOUT = RATOUT - Ksb(l)
        END IF 
C43-----PRINT STREAMFLOWS AND RATES FOR EACH REACH TO STREAM LIST
        IF (Kkstp==NSTP(KKPER)) THEN
          WRITE (912, 9010) il, ir, ic, ISTRMdfw(4,l), ISTRMdfw(5,l),  
     +                   Q(l), -1.0*Ksb(l), Stor1, STRMdfw(9,l),                         
     +                   Hnow, SNGL(depth), STRMdfw(5,l), P,         
     +                   STRMdfw(6,l), KSTRM(IA(l+1)-2),SNGL(hsb), kkper                    
        END IF
!MAS 13/1/15 ADDED PRINT STATEMENT FOR TIMING CALIBRATION FOR ELIANA
           IF (kkstp == NSTP(KKPER)) THEN
             IF(KKPER>1) THEN
               pern = KKPER-1
               IF(istsg==tseg(pern)) THEN
                 IF(nreach==trch(pern)) THEN
               OPEN(920,FILE='DFW_ELIANA.txt')
               WRITE(920,9015)istsg,nreach,Q(l),STRMdfw(9,l),SNGL(depth)
                 END IF
               END IF
              END IF
           END IF 
      END DO
C49-----MOVE RATES, VOLUMES, AND LABELS INTO ARRAYS FOR PRINTING
      VBVL(3, MSUM) = RATIN
      VBVL(4, MSUM) = RATOUT
      VBVL(2, MSUM) = VBVL(2, MSUM) + RATOUT*DELT
      VBVL(1, MSUM) = VBVL(1, MSUM) + RATIN*DELT
      VBNM(MSUM) = text
      MSUM = MSUM + 1
C45-----RECORD STREAM GAGING STATION DATA
!      IF ( Iunitgage.GT.0 ) THEN
!	  rtime = TOTIM - DELT + irt*deltinc  !this needs to be fixed
!        CALL SGWF2GAG7DF(rtime, BUFF, ibd, Q)
!      END IF
C
 9009 FORMAT (1X, ///1X, A, '   PERIOD ', I6, '   STEP ', I8, //,
     +        ' LAYER ROW COL. STREAM  RCH.   FLOW INTO   ', 
     +        ' STREAM      DELTA     FLOW OUT OF' 
     +        '    STREAM    STREAM       STREAM      WETTED    ', 
     +        '  STREAMBED   STREAM    GW HEAD', /16X, 
     
     +        'SEG.NO.  NO.   STRM. RCH.   LEAKAGE   STORAGE     ', 
     +        'STRM. RCH.   ', 
     +        '  HEAD      DEPTH        WIDTH        PERI      ', 
     +        'CONDCTNCE    CONDUCT   ', /)
 9010 FORMAT (1X, I3, I5, I5, 2I6, 3X, 4(2X, 1E10.4), 2X, F10.5,
     +        5(2X, E10.4), 1X, F9.5, 1X, 3I)
 9015 FORMAT (1X, 2I6, 2(2X,1E10.4), 2X, E10.4)
      RETURN      
      END SUBROUTINE GWF2DFW1BD
!-----------------------------------------------------------------
      DOUBLE PRECISION FUNCTION dfwB(ij, Hx, iss)
!     Return value of DA flow equation
		 USE GWFDFWMODULE, ONLY:STRMdfw, KSTRM, ISTRMdfw, BRHS, KSB, HCOEF,
     +                      SEGDFW, ISEGDFW, IOTSG
      USE GLOBAL,       ONLY:HNEW, IOUT
      USE GWFNWTMODULE, ONLY:IA, JA
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
!     ARGUMENTS
!     -----------------------------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: CONFUN
      INTEGER ij, iss
		 DOUBLE PRECISION Hx, fflow
		 INTRINSIC SQRT, ABS
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
		 INTEGER jj, il, ir, ic, jseg, ireach
		 DOUBLE PRECISION Hsb, stor1, fkstrm, bc, coef
!     -----------------------------------------------------------------   
      dfwB = 0.0D0
      stor1 = 0.0D0
	fflow = 0.0D0
	coef = 0.0D0
!      OPEN(913,file='DFW_debug')
!  Add all stream terms first; Kstrm*Hconn  - Kstrm*Hx
      Do jj = IA(ij)+1,IA(ij+1)-2 
        fflow =  KSTRM(jj)*sqrt(abs(STRMdfw(15,JA(jj))-Hx))
	  IF ((STRMdfw(15,JA(jj))-Hx).LT.0) THEN
		 coef = -1.0
	  ELSE 
		 coef = 1.0
        END IF
        dfwB = dfwB + coef*fflow
!		 WRITE(913,*)ij,JA(jj), KSTRM(JJ),HX, dfwB
      END DO
! ADD GW
	dfwB = dfwB + Ksb(ij)
!US BC
	IF (ISTRMdfw(5,ij)==1)dfwB = dfwB + SEGDFW(2,ISTRMdfw(4, ij))
!MAS 8/23/10 ISEG(4, jseg is the last REACH in the segment, this was counting wrong
         jseg = ISTRMdfw(4, ij)
	   ireach = ISTRMdfw(5, ij)  
! Sets DS BC to streambed elevation
      IF (ireach.EQ.ISEGDFW(4, jseg).and.IOTSG(jseg).EQ.0)  THEN
	  BC = STRMdfw(3,ij)	  
	  IF ((BC-Hx).LT.0) THEN
		 coef = -1.0
		 fKstrm = CONFUN(ij, ij, ij, Hx, 0)  
	  ELSE 
		 coef = 1.0
		 fKstrm = CONFUN(ij, ij, ij, BC, 0)  
       END IF
       fflow = fKstrm*sqrt(abs(BC-Hx))
	  dfwB = dfwB + coef*fflow
      END IF
! Sets DS BC as flowin = flowout in last reach
!      IF (ireach.EQ.ISEGDFW(4, jseg).and.IOTSG(jseg).EQ.0)  THEN
!	  dfwB = dfwB -coef*fflow
!	END IF 
!
      END FUNCTION dfwB
!-----------------------------------------------------------------------
      SUBROUTINE STORFM(Hr, Hc, Hnow, Hold, l)
C     ******************************************************************
C     COMPUTE THE STORAGE TERM
C     VERSION  1.0.00: 2010
C     ******************************************************************
      USE GWFDFWMODULE, ONLY: STRMdfw,ISTRMdfw,ISEGDFW,NSTRMdfw,QSTAGE,
     +                        BRHS, SNDFW, SODFW, QTOL, HEPS
      USE GWFBASMODULE, ONLY:DELT
		 USE GWFNWTMODULE, ONLY:BB,A,IA
		 USE GLOBAL, ONLY:IOUT
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
!     -----------------------------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: LOOKUP, Sat_thickdfw, Dhorz
		 INTRINSIC ABS
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
!      INTEGER KPER
!      DIMENSION 
C     ------------------------------------------------------------------
C     LOCAL VARIABLES
      DOUBLE PRECISION Bk, Bkp1, Hnow, rchx, Term1, Term2, Htop, Hold
	DOUBLE PRECISION BT, BBT, dSdh, dmax, Hnow1, Hold1, Hr, Hc
	DOUBLE PRECISION Bkd,Bkp1d,dBkp,Term3, Sn, So
	INTEGER ireach, jseg, nstrpts, IDIM, l
C     ------------------------------------------------------------------
!      INTEGER 
C     ------------------------------------------------------------------
!     Call interpolation tables to get width
      jseg = ISTRMdfw(4, l)         
	nstrpts = ISEGDFW(2, jseg)      
      IDIM = nstrpts+1     !B goes from nstrpts+1 to 2*nstrpts 
      Hold1 = Hold - STRMdfw(3,l)
      Hnow1 = Hnow - STRMdfw(3,l)
      Bk =  LOOKUP (Hold1, jseg, IDIM, nstrpts)
	Bkp1 = LOOKUP (Hnow1, jseg, IDIM, nstrpts) 
	Bkp1d = LOOKUP (Hnow1+HEPS, jseg, IDIM, nstrpts) 
	dBkp = (Bkp1d-Bkp1)/HEPS
      rchx = STRMdfw(1, l)    
!
	BT = STRMdfw(3, l)
	Htop = QSTAGE(nstrpts,jseg) + STRMdfw(3, l)
	jseg = ISTRMdfw(4, l)
	nstrpts = ISEGDFW(2, jseg)  
      Sn = Sat_thickdfw(Hnow, Htop, BT, l)
      So = Sodfw(l)
      dSdh = DHORZ(Hnow, Htop, BT)
	dmax = Htop - STRMdfw(3,l)
!		 
	Term1 = dmax*rchx/(DELT)
      Term2 = dBkp*Sn+dSdH*Bkp1
      Term3 = Bkp1*Sn-So*Bk
      Hc = Term1*Term2
      Hr = Term1*(Term2*Hnow-Term3)
	A(IA(l)) = A(IA(l)) - Hc
	BB(l) = BB(l) + Term1*Term3 - Hc*Hnow
	if(l.eq.1)then
!	write(iout,111)Hnow1,Bk,Bkp1,dBkp,Sn,So,term1*term2,term1*term3
	end if
111   format(9e20.10)
      END
!
!-----------------------------------------------------------------------
      SUBROUTINE STORBD(Hr, Hc, Hnow, Hold, l)
C     ******************************************************************
C     COMPUTE THE STORAGE TERM
C     VERSION  1.0.00: 2010
C     ******************************************************************
      USE GWFDFWMODULE, ONLY: STRMdfw,ISTRMdfw,ISEGDFW,NSTRMdfw,QSTAGE,
     +                        SNDFW, SODFW, QTOL
      USE GWFBASMODULE, ONLY:DELT
	USE GLOBAL, ONLY:IOUT
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
!     -----------------------------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: LOOKUP
	INTRINSIC ABS
C     ------------------------------------------------------------------
C     LOCAL VARIABLES
      DOUBLE PRECISION Bk, Bkp1, Hnow, rchx, Term1, Term2, Htop, Hold
		 DOUBLE PRECISION dmax, Hnow1, Hold1, Hr, Hc
		 DOUBLE PRECISION Sn, So
		 INTEGER ireach, jseg, nstrpts, IDIM, l
C     ------------------------------------------------------------------
      jseg = ISTRMdfw(4, l)         
	nstrpts = ISEGDFW(2, jseg)      
	IDIM = nstrpts+1     !B goes from nstrpts+1 to 2*nstrpts 
      Hold1 = Hold - STRMdfw(3,l)
      Hnow1 = Hnow - STRMdfw(3,l)
      Bk =  LOOKUP (Hold1, jseg, IDIM, nstrpts)
	Bkp1 = LOOKUP (Hnow1, jseg, IDIM, nstrpts) 
      rchx = STRMdfw(1, l)
      Htop = QSTAGE(nstrpts,jseg) + STRMdfw(3, l)    
!
	jseg = ISTRMdfw(4, l)
	nstrpts = ISEGDFW(2, jseg)  
      Sn = Sndfw(l)
      So = Sodfw(l)
      dmax = Htop - STRMdfw(3,l)
!		 
      Term1 = dmax*rchx/(DELT)
      Hc = Term1*Bkp1*Sn
      Hr = Term1*So*Bk
      END
!-----------------------------------------------------------------------
!     FUNCTION CONFUN, THE HORIZONTAL CONDUCTANCE  
      DOUBLE PRECISION FUNCTION CONFUN(ijup,ijnow,ijcon,Hnow,kflag)
	USE GWFDFWMODULE,  ONLY: ISTRMdfw, ISEGDFW, STRMdfw, QSTAGE,  
     +                         closezero, heps, CONST, sndfw
      USE GWFNWTMODULE,  ONLY: JA
	USE GLOBAL,        ONLY: iout
      IMPLICIT NONE
!     -----------------------------------------------------------------
!     ARGUMENTS
!     -----------------------------------------------------------------
      INTEGER  ijup, ijnow, ijcon, kflag
	DOUBLE PRECISION Hnow
!     -----------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: LOOKUP, Sat_thickdfw
      INTRINSIC SQRT, ABS
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      INTEGER nstrpts, jseg, IDIM
      DOUBLE PRECISION Area, Fn,A53,P,rchx,dmax,hc,Hmin, Sat,
     +		                   Htop, BT, P23
!     -----------------------------------------------------------------     
	jseg = ISTRMdfw(4, ijup)             !for nstrpts
	nstrpts = ISEGDFW(2, jseg)
	dmax = QSTAGE(nstrpts,jseg)          !max H for seg 
	IF (kflag == 0) hc = dmax*Sndfw(ijup)
	IF (kflag == 1) THEN
	  HTop = QSTAGE(nstrpts,jseg) + STRMdfw(3, ijup)
        BT = STRMdfw(3, ijup)
	  Sat = Sat_Thickdfw(Hnow, Htop, BT,ijup)
	  hc = dmax*Sat
	END IF
! call interpolation tables based on upstream head
! Area
      IDIM = 3*nstrpts+1        
      Area =  LOOKUP (hc, jseg, IDIM, nstrpts) 
! wetted perimeter
      IDIM = 4*nstrpts+1        
      P =  LOOKUP (hc, jseg, IDIM, nstrpts)
		 P23 = P**(2.0/3.0)
! manning n
      IDIM = 2*nstrpts+1      
      Fn =  LOOKUP (hc, jseg, IDIM, nstrpts) 
      A53 = Area**(5.0/3.0)
      IF (A53/P23.GT.CLOSEZERO) THEN
        rchx = (STRMdfw(1, ijnow)+STRMdfw(1, ijcon))*0.5
        CONFUN = (CONST/Fn)*(A53/P23)*(sqrt(rchx))/rchx
      ELSE
        CONFUN = 0.0
      END IF		 
!      write(iout,*) ijnow, hc, Area, P, confun
! 2000 format (3I,4(1x, E20.10))
      END FUNCTION CONFUN
!-----------------------------------------------------------------------
!     FUNCTION SBCON; Kvert LEAKAGE TO THE AQUIFER IN UNITS L3/T
      DOUBLE PRECISION FUNCTION SBCON(ij,Hchan,Hsb, Hup, iss)  
	USE GWFDFWMODULE,  ONLY: ISTRMdfw, ISEGDFW, STRMdfw, QSTAGE, 
     +                         SNDFW, SODFW, closezero, KSB, NSTRMdfw,
     +		 		 		 		 		 		  WETTIME
      USE GWFNWTMODULE,  ONLY: JA
	USE GWFBASMODULE, ONLY: DELT
	USE GLOBAL,        ONLY: iout
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     ARGUMENTS
!     ------------------------------------------------------------------
      INTEGER  ij, iss
	DOUBLE PRECISION Hsb
!     -----------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: LOOKUP, SAT_THICKDFW, philip, func9
	INTRINSIC SQRT, ABS
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      INTEGER jjseg, nstrpts, jseg, IDIM, il, ir, ic
      DOUBLE PRECISION Area, P, rchx, Hup, grad, Sat, a, b, seep,sp
      DOUBLE PRECISION Hchan, Htop, BT, Sn, hc, Kmax, Ks, sc,ss
!     -----------------------------------------------------------------
      jseg = ISTRMdfw(4, ij)             
	nstrpts = ISEGDFW(2, jseg)       
! Max stream stage
      HTop = QSTAGE(nstrpts,jseg) + STRMdfw(3, ij)
      BT = STRMdfw(3, ij)
	Sat = Sat_Thickdfw(Hup, Htop, BT,ij)
	hc = (Htop-STRMdfw(3, ij))*Sat
! Calc wetted perimeter for Hc    
      IDIM = 4*nstrpts+1      
      P =  LOOKUP (hc, jseg, IDIM, nstrpts)		   
	rchx = STRMdfw(1, ij)
	Ks = STRMdfw(6,ij)
	seep = Ks*DELT
	grad = (Hsb-Hchan)/STRMdfw(8,ij)		 
	IF (grad .LT. -1.0D0.AND.iss==0.AND.Ks.GT.0.0) THEN
	  grad = -1.0D0
	  ss = FUNC9(ij)
	  sc = sqrt(2.0*(STRMdfw(32,ij)-STRMdfw(33,ij))*(ss+Ks*hc))                                   
        a = WETTIME(ij) - delt                                  
        IF( a.LT.0.0 ) a = 0.0                            
        b = WETTIME(ij)                           
        seep = philip(b,a,seep,sc)                       
        SBCON = -1.0*(seep/delt)*rchx*P                       
	ELSE
	  IF (grad .GT. 1.0) grad = 1.0D0   
	  SBCON = Ks*rchx*P*grad 
	END IF
      END Function SBCON
!-----------------------------------------------------------------------
C     Philip's equation for vertical infiltration,from gwf2sfr7_philip.f      
      DOUBLE PRECISION FUNCTION philip(t,told,seep,sc)   
      DOUBLE PRECISION t,told,seep,sc,flobot             
      flobot = (sc*t**0.5 - sc*told**0.5 + seep)         
      IF (flobot.LT.0.0) flobot=0.0                      
      philip = flobot                                    
      END FUNCTION  
!-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FUNC9(ij)  
!this is the integrand from Brooks Corey = sorptivity
!from Wes Henson
!Rp=relative permeability                                
	USE GWFDFWMODULE, ONLY: STRMdfw                   
      DOUBLE PRECISION:: a, b, term1, term2, RP1,RP2,RP3,Ks
      RP1 =STRMdfw(34,ij)            !theta r        
      RP2 = STRMdfw(32,ij)		   !theta s  
      RP3 = 3.5                      !epsilon              
      Ks = STRMdfw(6,ij)		 	   !Ks                       
	a = STRMdfw(32,ij)             !theta sat
	b = STRMdfw(33,ij)             !theta init   
!THO term                                                         
	term1=((((RP1-a)/(RP1-RP2))**RP3*(a-RP1))/(RP3+1))    
!THI term                                                         
	term2=((((RP1-b)/(RP1-RP2))**RP3*(b-RP1))/(RP3+1))    
	func9=Ks*(term1-term2)                                      
      END FUNCTION FUNC9    
!-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION SAT_THICKDFW(Hup,Ttop,Bbot,ij)
! RETURNS SATURATED THICKNESS OF CELL BASED ON SMOOTH FUNCTION
      USE GWFDFWMODULE
      USE GLOBAL, ONLY: IOUT
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      INTEGER ij
      DOUBLE PRECISION hup, bbot, zero, ttop, factor, x, EPS, ACOF 
      DOUBLE PRECISION  Y
!     -----------------------------------------------------------------
      zero = 0.0D0
      SAT_THICKDFW = 1.0D0
!    STRAIGHT LINE WITH PARABOLIC SMOOTHING
      EPS = Thickfactdfw
      ACOF = 1.0 / (1.0 - EPS)
      x = (Hup-bbot)/(TTOP-BBOT)
      IF ( x.LT.1.0e-9 )x = 1.0e-9
      IF ( x.LT.0.0 ) THEN
        Y = 0.0
      ELSEIF(X.LT.EPS)THEN
        Y = ACOF *0.5/EPS * X**2
      ELSEIF(X.LT.1.0-EPS)THEN
        Y = ACOF * X + (1.0-ACOF)*0.5
      ELSEIF(X.LT.1.0)THEN
        X = 1.0 - X
        Y = ACOF *0.5/EPS * X**2
        Y = 1.0-Y
      ELSE
        Y = 1.0
      ENDIF
      factor = Y
      SAT_THICKDFW = factor
      END FUNCTION SAT_THICKDFW
!-----------------------------------------------------------------------
!     Updates saturation for latest iteration.
      SUBROUTINE Sndfw_update()
      USE GWFDFWMODULE, ONLY: NSTRMDFW,ISTRMDFW,ISEGDFW,QSTAGE,STRMDFW,
     +                        SNDFW
      USE GLOBAL,      ONLY:Iout
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: Sat_thickdfw
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      DOUBLE PRECISION HH, TP, BT, BBT
      INTEGER ij, ireach, jseg, nstrpts
!     -----------------------------------------------------------------   
      DO ij = 1, NSTRMdfw
        jseg = ISTRMdfw(4, ij)
	  nstrpts = ISEGDFW(2, jseg)  
        HH = STRMdfw(15, ij)
        TP = QSTAGE(nstrpts,jseg) + STRMdfw(3,ij)
        BT = STRMdfw(3, ij)
        Sndfw(ij) = Sat_thickdfw(HH, TP, BT, ij)
      END DO
      END SUBROUTINE
!-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DHORZ(Hup, Ttop, Bbot)
! RETURNS DERIV OF HORIZ K (dSn/dh) BASED ON SMOOOOOTH FUNCTION, MODIFIED 
! FROM THE NWT SUBROUTINE DHORIZ
      USE GWFDFWMODULE
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     ARGUMENTS
!     -----------------------------------------------------------------
      DOUBLE PRECISION Hup, Ttop, Bbot, Bbotm
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      DOUBLE PRECISION factor, x, EPS, ACOF, Y
!     -----------------------------------------------------------------
C-------STRAIGHT LINE WITH PARABOLIC SMOOTHING
      EPS = Thickfactdfw
      ACOF = 1.0 / (1.0 - EPS)
      x = (Hup-bbot)/(TTOP-BBOT)
	IF ( x.LT.1.0e-9 )x = 1.0e-9 
      IF ( x.LT.0.0 ) THEN
        Y = 0.0
      ELSEIF(X.LT.EPS)THEN
        Y = ACOF * X / (EPS*(Ttop-Bbot))
      ELSEIF(X.LT.1.0D0-EPS)THEN
        Y = ACOF /(Ttop-Bbot)
      ELSEIF(X.LT.1.0D0)THEN
        X = 1.0 - X
        Y = - ACOF * x / (EPS * (Ttop - Bbot))
        Y = -Y
      ELSE
        Y = 0.0
      ENDIF
      factor = Y
      DHORZ = factor     
      END FUNCTION DHORZ
!
!-----------------------------------------------------------------------
!     FILL JACOBIAN
      SUBROUTINE Jacobyan(ISS)            
      USE GWFNWTMODULE, ONLY: JA, IA, A, BB, Hchange
	USE GWFDFWMODULE
	USE GLOBAL, ONLY: Iout
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
      INTEGER ISS
      DOUBLE PRECISION, EXTERNAL ::SBCON,BBHeps, BbOffdiag, dfwB
	INTRINSIC SQRT, ABS
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     ----------------------------------------------------------------- 
      DOUBLE PRECISION Hnow, Hnower1, Hnower2, rms, ferr, tolf, Adiag
	DOUBLE PRECISION Hsb, Hold, Hgw1, Hgw2, zero, Hcon, temp, Hup
      DOUBLE PRECISION Hoffdi1, Hoffdi2, BC, stor1, test, sum, dHcoef 
	DOUBLE PRECISION Hcoefheps, Brhsheps, dBrhs, fKsb1, fKsb2
	DOUBLE PRECISION Hmeps, Hpeps
      INTEGER ij, l, irchfld, isegfld, jj, ijup, jseg
!     -----------------------------------------------------------------
      tolf = 1.0D-12
	rms1dfw = 0.0D0
	rms = 0.0D0
	dHcoef = 0.0D0   
	Hcoefheps = 0.0D0  
	Brhsheps = 0.0D0  
	dBrhs = 0.0D0  
	fKsb1 = 0.0D0 
	fKsb2 = 0.0D0 
	stor1 = 0.0D0  
 ! For each row (reach)
      DO ij = 1, NSTRMdfw  
        sum  = 0.0D0 
	  ferr = func(ij)          
	  rms = rms + ferr**2
	  Hsb = STRMdfw(19, ij)
	  Hold = STRMdfw(31, ij)
! DIAGONAL
	  Hnow = STRMdfw(15, ij)     		 		   
        Hnower1 = Hnow + HEPS
	  Hpeps = BbHeps(Hnow, Hnower1, ij, ISS)
	  Hnower2 = Hnow
	  Hmeps = dfwB(ij, Hnow, iss)
		   !!
	  l = IA(ij)
        test = (Hpeps - Hmeps)/(HEPS)
        A(l) = A(l) + test
	  sum = sum + A(l)   !for linear solver tolerance
	  BB(ij) = BB(ij) - Hmeps + test*Hnow   !right hand side
!		 write(iout,9001)ij,Hnow, Func(ij), Hpeps, Hmeps
!     +      BbB(ij), A(l)
!      write(iout, *) 'diag', ij, func(ij), A(l), Hpeps

! GROUNDWATER OFFDIAGONAL		     
	  Hgw1 = Hsb + HEPS
        If (Hnow .GT. Hgw1) THEN
		Hup = Hnow
	  ELSE
		Hup = Hgw1
	  END IF 
!	  Hgw2 = Hsb
	  fKsb1 = SBCON(ij, Hnow, Hgw1, Hup, iss)
!	  fKsb2 = SBCON(ij, Hnow, Hgw2, Hnow, iss)  
	  dgwQ(ij) = (fKsb1 - Ksb(ij))/(HEPS)  
	  l = IA(ij+1)-1          
	  A(l) = A(l) + dgwQ(ij) 
	  sum = sum + dgwQ(ij)   !for linear solver tolerance
	  BB(ij) = BB(ij) + dgwQ(ij)*Hsb   !right hand side
!		 write(iout,9003) ij, JA(l), fKsb, Ksb(ij), A(l)
!      write(iout, *) 'GW', ij, dqdhstr(ij),dgwQ(ij)
!		    
! OFFDIAGONALS (except the GW):
        DO l = IA(ij)+1,IA(ij+1)-2
		 Hnow = STRMdfw(15, ij)
		 Hoffdi1 = STRMdfw(15, JA(l)) + HEPS
		 Hpeps = BbOffdiag(Hnow, Hoffdi1, ij, l, iss)
		 Hoffdi2 = STRMdfw(15, JA(l))
!		    Hmeps = BbOffdiag(Hnow, Hoffdi2, ij, l, iss)
		 Test = (Hpeps - Hmeps)/(HEPS)                !!
           A(l) = A(l) + Test
	     sum = sum + Test   !for linear solver tolerance
		 BB(ij) = BB(ij) + Test*STRMdfw(15, JA(l))   !right hand side
!		   write(iout,9002)ij,JA(l),Hnow,Hcon,Func(ij),bbheps,fkstrm,
!     +                 Kstrm(ij)
!      write(iout, *) 'offdiag', JA(l), A(l), BB(ij)
        End Do
! Check for row values bad for linear solver
        Adiag = ABS(A(IA(ij)))
        IF ( Adiag.LT.tolf ) THEN
           A(IA(ij)) = 1.0D0
           BB(ij) = BB(ij) + 1.0D0*STRMdfw(3, ij)
        ELSE IF ( abs(sum).LT.tolf ) THEN
 !           A(IA(ij)) = A(IA(ij)) + tolf
 !           BB(ij) = BB(ij) + tolf*STRMdfw(3, ij)
           A(IA(ij)) = 1.0D0
           BB(ij) = BB(ij) + 1.0D0*STRMdfw(3, ij)
        END IF
        IF ( abs(ferr).GT.abs(Ffluxdfw) ) THEN
           ffluxdfw = ferr
           irchfld = ISTRMdfw(5, ij)
           isegfld = ISTRMdfw(4, ij)
        END IF
!Initial value for solver	  
        Hchange(ij) = STRMdfw(15, ij)   	  
      End Do 
      rms1dfw = rms**0.5
      rmsdfw = rms1dfw
  111 FORMAT ('Reach    Hnow                Func              BBheps '              
     + '                fKstrm               Kstrm                    '         
     + 'BbB                   A')
 9001 FORMAT (1X,I3,3X, E10.4, 6(3X,E20.10))
  112 FORMAT ('Reach  Conxn   Hnow    Hcon              Func          '              
     +'     BBheps                     fKstrm              Kstrm      ')         
 9002 FORMAT (1X,I3,2X,I3, 2(2X, E10.4), 4(2X,E20.10))
  113 FORMAT ('Reach  Conxn   fKsb               Ksb              A')                       
 9003 FORMAT (1X,I3,2X,I3, 2X, E10.4, 2(3X,E20.10))
      END SUBROUTINE Jacobyan
!-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION BbHeps(Hnow, Hnower, ij, ISS)
! RETURNS BBHEPS FOR THE NUMERICAL DERIVATIVE 
      USE GWFDFWMODULE
	USE GWFNWTMODULE, ONLY: JA, IA, A, BB, Hchange
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ARGUMENTS
!     -----------------------------------------------------------------
      INTEGER ISS
      DOUBLE PRECISION, EXTERNAL ::dfwB,CONFUN,SBCON
      INTRINSIC SQRT, ABS
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     ----------------------------------------------------------------- 
      DOUBLE PRECISION Hnow, Hnower, deriv, rms, ferr, tolf, Adiag, sum
      DOUBLE PRECISION Bboffdi, Hsb, Hold, Hgw2, zero, Hcon, temp, Hup
      DOUBLE PRECISION Hoffdi, Hoffdi2, BC, test, fflow
	DOUBLE PRECISION Hcoefheps, Brhsheps, fKsb, fKstrm, coef
      INTEGER ij, l, irchfld, isegfld, jj, ijup, jseg, kflag, ireach
!     -----------------------------------------------------------------
	BbHeps = 0.0D0  
	Hsb = STRMdfw(19, ij)
	Hold = STRMdfw(31, ij)  
	Hcoefheps = 0.0D0  
	Brhsheps = 0.0D0   
	fKstrm = 0.0D0 
	kflag = 1 
	Do jj = IA(ij)+1,IA(ij+1)-2
	  Hcon = STRMdfw(15,JA(jj))  		 
	  IF ((Hcon-Hnower).LT.0) THEN
		 coef = -1.0           
		 ijup = ij
		 fKstrm = CONFUN(ijup, ij, JA(jj), Hnower, kflag)
		 fflow =  fKSTRM*sqrt(abs(Hcon-Hnower))
		 BbHeps = BbHeps + coef*fflow
	  ELSE 
		 coef = 1.0
		 ijup = JA(jj)
		 fflow = KSTRM(jj)*sqrt(abs(Hcon-Hnower))
		 BbHeps = BbHeps + coef*fflow
        END IF     
      END DO

	If (Hnower .GT. Hsb) THEN
		Hup = Hnower
	ELSE
		Hup = Hsb
	END IF 
	fKsb = SBCON(ij, Hnower, Hsb, Hup, iss)  !GW
	BbHeps = BbHeps + fKsb
!Change in leakage with respect to SW head- only for NWT offdiag
	dqdhstr(ij) = (fKsb - Ksb(ij))/(HEPS) 
      IF (ISTRMdfw(5,ij)==1) 
     +    BbHeps = BbHeps + SEGDFW(2,ISTRMdfw(4, ij))  !Inflow
	jseg = ISTRMdfw(4, ij) 
	ireach = ISTRMdfw(5, ij)  !MAS added 8/23/10
! Sets DS BC to streambed elevation
	IF (ireach.EQ.ISEGDFW(4, jseg).and.IOTSG(jseg).EQ.0)  THEN  
	  ijup = ij
	  BC = STRMdfw(3,ij)
	  IF ((BC-Hnower).LT.0) THEN
		  coef = -1.0
		  fKstrm = CONFUN(ijup, ij, ij, Hnower, kflag)
	  ELSE 
		  coef = 1.0
		  fKstrm = CONFUN(ijup, ij, ij, BC, kflag)
        END IF
	  fflow = fKstrm*sqrt(abs(BC-Hnower))
	  BbHeps = BbHeps + coef*fflow
	END IF
! Sets DS BC as flowin = flowout in last reach
!      IF (ireach.EQ.ISEGDFW(4, jseg).and.IOTSG(jseg).EQ.0)  THEN
!	  BbHeps = BbHeps -coef*fflow
!	END IF 
!		 write(iout,*) Brhs(ij), BrhsHeps, BbHeps
      END FUNCTION
!------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION Bboffdiag(Hnow, Hoffdi, ij, l, ISS)
! RETURNS BBHEPS FOR THE NUMERICAL DERIVATIVE FOR THE OFFDIAG
      USE GWFDFWMODULE
		 USE GWFNWTMODULE, ONLY: JA, IA, A, BB, Hchange
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ARGUMENTS
!     -----------------------------------------------------------------
      DOUBLE PRECISION, EXTERNAL ::dfwB,CONFUN,SBCON
		 INTEGER iss
		 INTRINSIC SQRT, ABS
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     ----------------------------------------------------------------- 
      DOUBLE PRECISION Hnow, Hnower, deriv, rms, ferr, tolf, Adiag, sum
	DOUBLE PRECISION Bboffdi, Hsb, Hold, Hgw2, zero, Hcon, temp, Hup
      DOUBLE PRECISION Hoffdi, Hoffdi2, BbHeps, BC, test, fflow
	DOUBLE PRECISION Hcoefheps, Brhsheps, fKsb, fKstrm,coef
	DOUBLE PRECISION Hr, Hc
      INTEGER ij, l, irchfld, isegfld, jj, ijup, jseg, kflag, ireach
!     -----------------------------------------------------------------
	BbHeps = 0.0D0  
	Hsb = STRMdfw(19, ij)
	Hold = STRMdfw(31, ij)  
	Hcoefheps = 0.0D0  
	Brhsheps = 0.0D0 
	fKstrm = 0.0D0 
	kflag = 1
	Do jj = IA(ij)+1,IA(ij+1)-2
	  IF ( l.EQ.jj ) THEN 
		 IF ((Hoffdi-Hnow).LT.0) THEN
		    coef = -1.0
			ijup = ij
!		    fKstrm = CONFUN(ijup,ij,JA(jj),Hnow,kflag)
		    fflow = KSTRM(jj)*sqrt(abs(Hoffdi-Hnow))
		    BbHeps = BbHeps + coef*fflow
		 ELSE 
		    coef = 1.0              
			ijup = JA(l)
		    fKstrm = CONFUN(ijup,ij,JA(jj),Hoffdi,kflag)
		 	fflow = fKSTRM*sqrt(abs(Hoffdi-Hnow))
		    BbHeps = BbHeps + coef*fflow
           END IF       
        ELSE
           Hoffdi2 = STRMdfw(15, JA(jj))
           IF ((Hoffdi2-Hnow).LT.0) THEN
		    coef = -1.0
			ijup = ij
 !             fKstrm = CONFUN(ijup, ij, JA(jj), Hnow, kflag)
		    fflow = KSTRM(jj)*sqrt(abs(Hoffdi2-Hnow))
		    BbHeps = BbHeps + coef*fflow
		 ELSE 
		    coef = 1.0
			ijup = JA(jj)
              fKstrm = CONFUN(ijup, ij, JA(jj), Hoffdi2, kflag)
     	        fflow = fKSTRM*sqrt(abs(Hoffdi2-Hnow))
		    BbHeps = BbHeps + coef*fflow
           END IF     
	  END IF
      END DO
	BbHeps = BbHeps + Ksb(ij)   !GW
	IF (ISTRMdfw(5,ij)==1) 
     +      BbHeps = BbHeps + SEGDFW(2,ISTRMdfw(4, ij))  !inflow                       
	jseg = ISTRMdfw(4, ij) 
	ireach = ISTRMdfw(5, ij)  !MAS added 8/23/10
!Sets DSBC to streambed elevation
	IF (ireach.EQ.ISEGDFW(4, jseg).and.IOTSG(jseg).EQ.0)  THEN  
	  ijup = ij
	  BC = STRMdfw(3,ij)
	  IF ((BC-Hnow).LT.0) THEN
	    coef = -1.0
	    fKstrm = CONFUN(ijup, ij, ij, Hnow, kflag)
	  ELSE 
	    coef = 1.0
	    fKstrm = CONFUN(ijup, ij, ij, BC, kflag)
        END IF
	  fflow = fKstrm*sqrt(abs(BC-Hnow))
		     BbHeps = BbHeps + coef*fflow
	END IF
! Sets DS BC as flowin = flowout in last reach 
!      IF (ireach.EQ.ISEGDFW(4, jseg).and.IOTSG(jseg).EQ.0)  THEN
!	  BbHeps = BbHeps -coef*fflow
!	END IF 
!
      BbOffdiag = BbHeps
      END FUNCTION
!------------------------------------------------------------------
C-------SUBROUTINE SGWF2DFW1RDSEG
      SUBROUTINE SGWF2DFW1RDSEG(Nlst, Lstbeg, In, Ischk, 
     +                          Nischk, Ichk, Kkper)
C     ******************************************************************
C     READ STREAM SEGMENT DATA 
C     VERSION  1.0.00: 2009
C     ******************************************************************
      USE GWFDFWMODULE, ONLY: NSS, MAXPTS, IDIVA, IOTSG, ISEGDFW,
     +                        SEGDFW,  QSTAGE 
      USE GLOBAL,       ONLY: IOUT
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER Ichk, In, Ischk, Lstbeg, Nischk, Nlst, Kkper
      DIMENSION Ischk(Nischk)
C     ------------------------------------------------------------------
C     LOCAL VARIABLES
C     ------------------------------------------------------------------
      INTEGER icalc, idum, ii, iqseg, iupseg, jj, jk, lstend, n, 
     +        noutseg, nseg, nstrpts
C     ------------------------------------------------------------------
C
C1------READ STREAM SEGMENT DATA.
      lstend = Lstbeg + Nlst - 1
      DO iqseg = Lstbeg, lstend
C
C2------ONLY READ FIRST 4 VARIABLES TO DETERMINE VALUE OF IUPSEG.
        READ (In, *) n, icalc, noutseg, iupseg
        IF ( n.GT.NSS .OR. n.LT.1 ) THEN
          WRITE (IOUT, 9001) n
 9001     FORMAT (1X, /1X, 'SEGMENT NUMBER (NSEG) OUT OF RANGE: ', I6)
          IF ( Ichk.NE.0 ) THEN
            WRITE (IOUT, 9002) iqseg - Lstbeg + 1
 9002       FORMAT (1X, 'READING ENTRY ', I6, ' OF ITEM 6A')
          ELSE
            WRITE (IOUT, 9003) iqseg - Lstbeg + 1
 9003       FORMAT (1X, 'READING ENTRY ', I6, ' OF ITEM 4A')
          END IF
          CALL USTOP(' ')
        END IF
C
C3------STORE DATA IN ACTIVE SEGMENT AREA 
        nseg = n
        Ischk(n) = Ischk(n) + 1
        BACKSPACE In
C
C4------READ DATA SET 4B FOR SEGMENTS THAT ARE NOT DIVERSIONS.
        IF ( iupseg.LE.0 ) THEN
            READ (In, *) idum, ISEGDFW(1, nseg), IOTSG(nseg), 
     +                   IDIVA(1, nseg), ISEGDFW(2, nseg), 
     +                   (SEGDFW(jj, nseg), jj=2, 5)
C
C5------READ DATA 4B FOR SEGMENTS THAT ARE DIVERSIONS FROM STREAMS.
        ELSE 
          READ (In, *) idum, ISEGDFW(1, nseg), IOTSG(nseg), 
     +                 (IDIVA(ii, nseg), ii=1, 2), ISEGDFW(2, nseg), 
     +                 (SEGDFW(jj, nseg), jj=2, 5)
        END IF
C
C
C9------READ DATA SET 4F FOR SEGMENT WHEN ICALC IS 4.
        IF ( icalc.EQ.4 ) THEN
          nstrpts = ISEGDFW(2, nseg)
          IF ( nstrpts.LT.2 ) THEN
            WRITE (IOUT, 9004) n
 9004       FORMAT (/1X, 'NUMBER OF POINTS USED TO RELATE ', 
     +              'STREAMFLOW WITH STREAM DEPTH AND WIDTH FOR ', 
     +              'SEGMENT ', I6, ' IS LESS THAN TWO'//1X, 
     +              'PROGRAM STOPPING')
            CALL USTOP(' ')
          ELSE IF ( nstrpts.GT.MAXPTS/5 ) THEN
            WRITE (IOUT, 9005) n, nstrpts
 9005       FORMAT (/1X, 'FOR SEGMENT ', I6, ' NUMBER OF POINTS', 
     +              'USED TO RELATE STREAMFLOW WITH DEPTH AND ', 
     +              'WIDTH IS ', I5//1X, 'WHICH IS MORE THAN ', 
     +              'MAXIMUM NUMBER OF 50 POINTS', //1X, 
     +              'PROGRAM STOPPING'//)
            CALL USTOP(' ')
          ELSE
            READ (In, *) (QSTAGE(jj, nseg), jj=1, nstrpts)             !h
            READ (In, *) (QSTAGE(jj, nseg), jj=nstrpts+1, 2*nstrpts)   !b
		  READ (In, *) (QSTAGE(jj, nseg), jj=2*nstrpts+1, 3*nstrpts) !n
		  READ (In, *) (QSTAGE(jj, nseg), jj=3*nstrpts+1, 4*nstrpts) !A
		  READ (In, *) (QSTAGE(jj, nseg), jj=4*nstrpts+1, 5*nstrpts) !P
          END IF
        END IF
      END DO
!
      END SUBROUTINE SGWF2DFW1RDSEG
!-----------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION LOOKUP (STAGE, jseg, IDIM, nstrpts)
C     *******************************************************************
C     COMPUTE WIDTH AREA MANNING HYDRAD GIVEN STAGE USING RATING TABLES.
C     LINEARLY INTERPOLATES BETWEEN TWO VALUES, with optional smoothing curve
C     between the first two values.
C     VERSION  0.0.01 (FINTERP Donated from the RGN files) 9/2/09
C     *******************************************************************     
      USE GWFDFWMODULE, ONLY: QSTAGE, STRMDFW 
      IMPLICIT NONE              
      DOUBLE PRECISION FOLD, TOLF2, s, thick, Qp
      DOUBLE PRECISION STAGE,check, ttop, x, aa, b, cof1,cof2,cof3
      INTEGER I, IFLG, IDIM, imax, nstrpts, jseg, ij  
!------------------------------------------------------------------------
      TOLF2=1.0E-14
      imax = nstrpts
      LOOKUP = 0.0D0
!optional smoothing curve between first two points
!      IF (STAGE.GT.0.0.AND.STAGE.LT.QSTAGE(2,jseg)) THEN
!	  ttop = QSTAGE(2, jseg)      ! second point in the table
!        x = STAGE 
!        cof1 = x**2.0D0
!        cof2 = -(2.0D0*x)/(ttop**3.0D0)
!        cof3 = 3.0D0/(ttop**2.0D0)
!        Qp = cof1*(cof2+cof3)
!        LOOKUP = Qp*QSTAGE(IDIM+1,jseg)
!        RETURN
!      END IF
!INPUT = H values = QSTAGE(1:nstrpts)
!IF H greater than Hmax, Set to Hmax
      IF (STAGE.GT.QSTAGE(imax,jseg))THEN
        LOOKUP =  QSTAGE(imax+IDIM-1,jseg) 
        RETURN
!IF H below Hmin, keep the width at 0
      ELSE IF (STAGE.LT.QSTAGE(1,jseg)) THEN
        LOOKUP =  QSTAGE(IDIM,jseg)  
	  RETURN
      END IF
      IFLG = 0
      I = 1
      DO WHILE ( IFLG.EQ.0 )
        FOLD=ABS(STAGE-QSTAGE(I,jseg))     
        IF (FOLD .LE. TOLF2) THEN  
          LOOKUP =QSTAGE(IDIM,jseg)
          IFLG = 1
        ELSEIF (STAGE.GT.QSTAGE(I,jseg) .AND. STAGE.LT.
     +          QSTAGE(I+1,jseg))THEN
          LOOKUP =((QSTAGE(IDIM+I,jseg)-QSTAGE(IDIM + I-1,jseg))/
     +         (QSTAGE(I+1,jseg)- QSTAGE(I,jseg)))*
     +         STAGE+QSTAGE(IDIM+ I,jseg)-((QSTAGE(IDIM+I,jseg)-
     +         QSTAGE(IDIM+ I-1,jseg))/(QSTAGE(I+1,jseg)-
     +         QSTAGE(I,jseg)))*QSTAGE(I+1,jseg)                 
          IFLG = 1
        END IF
        I = I + 1
        IF( I.GT.nstrpts-1 ) IFLG = 1 
      END DO
      LOOKUP = LOOKUP  
	IF ( LOOKUP.LT.1.0e-12 ) LOOKUP = 1.0e-12
!
      END FUNCTION LOOKUP
!-----------------------------------------------------------------------
C-------SUBROUTINE SGWF2DFW1PRSEG
      SUBROUTINE SGWF2DFW1PRSEG(Nlst, Lstbeg, Kkper, Iouts)
C     ******************************************************************
C     PRINT STREAM SEGMENT DATA 
C     VERSION  1.0.00:  2009
C     ******************************************************************
      USE GWFDFWMODULE, ONLY: IDIVA, IOTSG, ISEGDFW, SEGDFW,
     +                        QSTAGE
      USE GLOBAL,       ONLY: IOUT
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER Nlst, Lstbeg, Kkper, Iouts
C     ------------------------------------------------------------------
C     LOCAL VARIABLES
C     ------------------------------------------------------------------
      INTEGER i, icalc, iflg, ii, ipt, jj, lstend, nn, nseg, 
     +        nstrpts
C     ------------------------------------------------------------------
C
      lstend = Nlst + Lstbeg - 1
      WRITE (IOUT, 9001)
 9001 FORMAT (1X, //20X, 'INPUT DATA FOR EACH STREAM SEGMENT', /1X, 
     +        93('-')/)
C
C1------PRINT INPUT FLOW RATES FOR EACH STREAM SEGMENT.
      WRITE (IOUT, 9002)
 9002 FORMAT (1X, 'SEGMENT    SEG.     INFLOW   OVERLAND   ', 
     +        'STREAM    STREAM   ICALC  OUTFLOW  DIVERSION PRIORITY', 
     +        /4X, 'NO.    LENGTH     RATE     RUNOFF      ', 
     +        'ET       PPT.    METH.  TO SEG.  FROM SEG.    NO.'/)
      DO nseg = Lstbeg, lstend
        IF ( Lstbeg.EQ.1 ) THEN
          nn = nseg
        ELSE
          nn = ISEGDFW(3, nseg)
        END IF
        WRITE (IOUT, 9003) nn,(SEGDFW(ii, nseg),ii=1,5),ISEGDFW(1,nseg), 
     +                     IOTSG(nseg), (IDIVA(jj, nseg), jj=1, 2)
 9003   FORMAT (1X, I6, 1X, 1P5E10.3, 2X, I3, 3X, I6, 3X, I6, 4X, I5)
      END DO
C
C2------PRINT STREAMBED PROPERTIES AND STREAM DIMENSIONS.
!      IF ( Lstbeg.EQ.1 ) THEN
!          WRITE (IOUT, 9005)
!      ELSE
!        WRITE(IOUT,210)
!  210   FORMAT (1X,//9X,'STREAMBED PROPERTIES AND STREAM ',
!     1        'DIMENSIONS',//1X,'SEGMENT  BED HYD. COND. FACTOR',2X,
!     2        'BED THICKNESS     ELEV.-TOP OF BED     WIDTH OF ',
!     3        'STREAM     DEPTH OF STREAM    STREAM ROUGHNESS',/1X,
!     4        '   No.     UPPER     LOWER     UPPER     ',
!     5        'LOWER     UPPER     LOWER     UPPER     LOWER     ',
!     6        'UPPER     LOWER   CHANNEL      BANK'/)
!      END IF
! 9005 FORMAT (1X, //9X, 'STREAMBED PROPERTIES AND STREAM DIMENSIONS', //
!     +        ' SEGMENT     WIDTH OF STREAM', 5X,
!     +        'DEPTH OF STREAM    STREAM ROUGHNESS', /, 
!     +        '    No.     UPPER     LOWER     UPPER     ', 
!     +        'LOWER     CHANNEL      BANK'/)
C
C6------PRINT TABULATED VALUES FOR COMPUTING STREAM WIDTH AND DEPTH
C         FROM STREAMFLOW FOR SEGMENTS WITH ICALC=4.
      iflg = 0
      DO nseg = Lstbeg, lstend
        IF ( Lstbeg.EQ.1 ) THEN
          nn = nseg
        ELSE
          nn = ISEGDFW(3, nseg)
        END IF
        icalc = ISEGDFW(1, nseg)
        nstrpts = ISEGDFW(2, nseg)
        IF ( icalc.EQ.4 .AND. iflg.EQ.0 ) THEN
          WRITE (IOUT, 9031)
 9031     FORMAT (1X, /1X, 'STREAMFLOW RELATION WITH DEPTH ', 
     +            'AND WIDTH IS BASED ON TABULATED VALUES', //2X, 
     +            'SEGMENT NO.   DEPTH        WIDTH        MANNING ', 
     +            '     AREA          WET PERIM', /)
          iflg = 1
        END IF
        ipt = 1
        IF ( icalc.EQ.4 .AND. iflg.EQ.1 ) THEN
          DO WHILE ( ipt.LE.nstrpts )
            WRITE (IOUT, 9032) nn, QSTAGE(ipt, nseg), 
     +                         QSTAGE(nstrpts+ipt, nseg), 
     +                         QSTAGE(2*nstrpts+ipt, nseg),
     +                         QSTAGE(3*nstrpts+ipt, nseg),
     +                         QSTAGE(4*nstrpts+ipt, nseg)
 9032       FORMAT (5X, I6, 2X, 5(3X, 1PE10.4))
            ipt = ipt + 1
          END DO
        END IF
      END DO
C
      RETURN
      END SUBROUTINE SGWF2DFW1PRSEG
!-----------------------------------------------------------------------
C-------SUBROUTINE GWF2DFW1DA
      SUBROUTINE GWF2DFW1DA(IGRID)
      USE GWFDFWMODULE
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER IGRID
C     ------------------------------------------------------------------
      DEALLOCATE (GWFDFWDAT(IGRID)%NSS)
      DEALLOCATE (GWFDFWDAT(IGRID)%NSTRMdfw)
      DEALLOCATE (GWFDFWDAT(IGRID)%ISTCB1)
      DEALLOCATE (GWFDFWDAT(IGRID)%ISTCB2)
      DEALLOCATE (GWFDFWDAT(IGRID)%MAXPTS)
      DEALLOCATE (GWFDFWDAT(IGRID)%Numactivedfw) 
      DEALLOCATE (GWFDFWDAT(IGRID)%ITMP)
      DEALLOCATE (GWFDFWDAT(IGRID)%IRDFLG)
      DEALLOCATE (GWFDFWDAT(IGRID)%IPTFLG)
      DEALLOCATE (GWFDFWDAT(IGRID)%IRTFLG)
      DEALLOCATE (GWFDFWDAT(IGRID)%CONST)
      DEALLOCATE (GWFDFWDAT(IGRID)%IOTSG)
      DEALLOCATE (GWFDFWDAT(IGRID)%NSEGCK)
      DEALLOCATE (GWFDFWDAT(IGRID)%ISEGDFW)
      DEALLOCATE (GWFDFWDAT(IGRID)%IDIVA)
      DEALLOCATE (GWFDFWDAT(IGRID)%ISTRMdfw)
      DEALLOCATE (GWFDFWDAT(IGRID)%SGOTFL)
      DEALLOCATE (GWFDFWDAT(IGRID)%SEGDFW)
      DEALLOCATE (GWFDFWDAT(IGRID)%STRMdfw)
      DEALLOCATE (GWFDFWDAT(IGRID)%QSTAGE)
      DEALLOCATE (GWFDFWDAT(IGRID)%dqdhstr)
		 DEALLOCATE (GWFDFWDAT(IGRID)%WETTIME)
		 DEALLOCATE (GWFDFWDAT(IGRID)%SNDFW)
		 DEALLOCATE (GWFDFWDAT(IGRID)%SODFW)
		 DEALLOCATE (GWFDFWDAT(IGRID)%DGWQ)
		 DEALLOCATE (GWFDFWDAT(IGRID)%Q)
		 DEALLOCATE (GWFDFWDAT(IGRID)%LEAK)
		 DEALLOCATE (GWFDFWDAT(IGRID)%DFDH)
		 DEALLOCATE (GWFDFWDAT(IGRID)%THICKFACTDFW)
		 DEALLOCATE (GWFDFWDAT(IGRID)%KSTRM)
		 DEALLOCATE (GWFDFWDAT(IGRID)%Func)
		 DEALLOCATE (GWFDFWDAT(IGRID)%BRHS)
		 DEALLOCATE (GWFDFWDAT(IGRID)%HCOEF)
		 DEALLOCATE (GWFDFWDAT(IGRID)%RMS1dfw)
		 DEALLOCATE (GWFDFWDAT(IGRID)%RMS2dfw)
		 DEALLOCATE (GWFDFWDAT(IGRID)%RMSdfw)
		 DEALLOCATE (GWFDFWDAT(IGRID)%ICELLDFW)
		 DEALLOCATE (GWFDFWDAT(IGRID)%IGWINDEX)
		 DEALLOCATE (GWFDFWDAT(IGRID)%Fheaddfw)
		 DEALLOCATE (GWFDFWDAT(IGRID)%Fheadsavedfw)
		 DEALLOCATE (GWFDFWDAT(IGRID)%Ffluxdfw)
C
      END SUBROUTINE GWF2DFW1DA
C
C-------SUBROUTINE GWF2DFW1PNT
      SUBROUTINE SGWF2DFW1PNT(IGRID)
      USE GWFDFWMODULE
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER IGRID
C     ------------------------------------------------------------------
      NSS=>GWFDFWDAT(IGRID)%NSS
      NSTRMdfw=>GWFDFWDAT(IGRID)%NSTRMdfw
      ISTCB1=>GWFDFWDAT(IGRID)%ISTCB1
      ISTCB2=>GWFDFWDAT(IGRID)%ISTCB2
      MAXPTS=>GWFDFWDAT(IGRID)%MAXPTS
      Numactivedfw=>GWFDFWDAT(IGRID)%Numactivedfw  
      ITMP=>GWFDFWDAT(IGRID)%ITMP
      IRDFLG=>GWFDFWDAT(IGRID)%IRDFLG
      IPTFLG=>GWFDFWDAT(IGRID)%IPTFLG
      CONST=>GWFDFWDAT(IGRID)%CONST
      IRTFLG=>GWFDFWDAT(IGRID)%IRTFLG
      IOTSG=>GWFDFWDAT(IGRID)%IOTSG
      NSEGCK=>GWFDFWDAT(IGRID)%NSEGCK
      ISEGDFW=>GWFDFWDAT(IGRID)%ISEGDFW
      IDIVA=>GWFDFWDAT(IGRID)%IDIVA
      ISTRMdfw=>GWFDFWDAT(IGRID)%ISTRMdfw
      SGOTFL=>GWFDFWDAT(IGRID)%SGOTFL
      SEGDFW=>GWFDFWDAT(IGRID)%SEGDFW
      STRMdfw=>GWFDFWDAT(IGRID)%STRMdfw
      QSTAGE=>GWFDFWDAT(IGRID)%QSTAGE
      dqdhstr=>GWFDFWDAT(IGRID)%dqdhstr
		 WETTIME=>GWFDFWDAT(IGRID)%WETTIME
      SNDFW=>GWFDFWDAT(IGRID)%SNDFW
		 DGWQ=>GWFDFWDAT(IGRID)%DGWQ
		 Q=>GWFDFWDAT(IGRID)%Q
		 LEAK=>GWFDFWDAT(IGRID)%LEAK
      DFDH=>GWFDFWDAT(IGRID)%DFDH
		 SODFW=>GWFDFWDAT(IGRID)%SODFW
		 THICKFACTDFW=>GWFDFWDAT(IGRID)%THICKFACTDFW
		 KSTRM=>GWFDFWDAT(IGRID)%KSTRM
		 Func=>GWFDFWDAT(IGRID)%Func
		 BRHS=>GWFDFWDAT(IGRID)%BRHS
		 HCOEF=>GWFDFWDAT(IGRID)%HCOEF
		 RMS1dfw=>GWFDFWDAT(IGRID)%RMS1dfw
		 RMS2dfw=>GWFDFWDAT(IGRID)%RMS2dfw
		 RMSdfw=>GWFDFWDAT(IGRID)%RMSdfw
		 ICELLDFW=>GWFDFWDAT(IGRID)%ICELLDFW 
		 IGWINDEX=>GWFDFWDAT(IGRID)%IGWINDEX  
		 Fheaddfw=>GWFDFWDAT(IGRID)%Fheaddfw
		 Fheadsavedfw=>GWFDFWDAT(IGRID)%Fheadsavedfw
		 Ffluxdfw=>GWFDFWDAT(IGRID)%Ffluxdfw
C
      END SUBROUTINE SGWF2DFW1PNT
C
C-------SUBROUTINE SGWF2DFW1PSV
      SUBROUTINE SGWF2DFW1PSV(IGRID)
      USE GWFDFWMODULE
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER IGRID
C     ------------------------------------------------------------------
      GWFDFWDAT(IGRID)%NSS=>NSS
      GWFDFWDAT(IGRID)%NSTRMdfw=>NSTRMdfw
      GWFDFWDAT(IGRID)%ISTCB1=>ISTCB1
      GWFDFWDAT(IGRID)%ISTCB2=>ISTCB2
      GWFDFWDAT(IGRID)%MAXPTS=>MAXPTS
      GWFDFWDAT(IGRID)%Numactivedfw =>Numactivedfw
      GWFDFWDAT(IGRID)%ITMP=>ITMP
      GWFDFWDAT(IGRID)%IRDFLG=>IRDFLG
      GWFDFWDAT(IGRID)%IPTFLG=>IPTFLG
      GWFDFWDAT(IGRID)%CONST=>CONST
      GWFDFWDAT(IGRID)%IRTFLG=>IRTFLG
      GWFDFWDAT(IGRID)%IOTSG=>IOTSG
      GWFDFWDAT(IGRID)%NSEGCK=>NSEGCK
      GWFDFWDAT(IGRID)%ISEGDFW=>ISEGDFW
      GWFDFWDAT(IGRID)%IDIVA=>IDIVA
      GWFDFWDAT(IGRID)%ISTRMdfw=>ISTRMdfw
      GWFDFWDAT(IGRID)%SGOTFL=>SGOTFL
      GWFDFWDAT(IGRID)%SEGDFW=>SEGDFW
      GWFDFWDAT(IGRID)%STRMdfw=>STRMdfw
      GWFDFWDAT(IGRID)%QSTAGE=>QSTAGE
      GWFDFWDAT(IGRID)%dqdhstr=>dqdhstr
	GWFDFWDAT(IGRID)%WETTIME=>WETTIME
	GWFDFWDAT(IGRID)%SNDFW=>SNDFW
	GWFDFWDAT(IGRID)%SODFW=>SODFW
      GWFDFWDAT(IGRID)%DGWQ=>DGWQ
	GWFDFWDAT(IGRID)%Q=>Q
	GWFDFWDAT(IGRID)%LEAK=>LEAK
	GWFDFWDAT(IGRID)%DFDH=>DFDH
	GWFDFWDAT(IGRID)%THICKFACTDFW=>THICKFACTDFW
	GWFDFWDAT(IGRID)%KSTRM=>KSTRM
	GWFDFWDAT(IGRID)%Func=>Func
	GWFDFWDAT(IGRID)%BRHS=>BRHS
      GWFDFWDAT(IGRID)%HCOEF=>HCOEF
	GWFDFWDAT(IGRID)%RMS1dfw=>RMS1dfw
	GWFDFWDAT(IGRID)%RMS2dfw=>RMS2dfw
	GWFDFWDAT(IGRID)%RMSdfw=>RMSdfw
	GWFDFWDAT(IGRID)%ICELLDFW=>ICELLDFW
	GWFDFWDAT(IGRID)%IGWINDEX=>IGWINDEX
	GWFDFWDAT(IGRID)%Fheaddfw=>Fheaddfw
	GWFDFWDAT(IGRID)%Fheadsavedfw=>Fheadsavedfw
	GWFDFWDAT(IGRID)%Ffluxdfw=>Ffluxdfw
C
      END SUBROUTINE SGWF2DFW1PSV