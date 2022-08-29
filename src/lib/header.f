       subroutine header(mpp,ipix,jpix)
C***********************************************************************
       character*1 cpl(4,256)
       character*4 byt4
       character*1 byta,bytb,bytc,bytd
       character*14 frmtstr
       character*54 headmsw
C-----------------rainbow color palette---------------------------------
C-----red
       do 110 i=1,142
       irgb=float(128*(29-i))/28.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(3,i)=CHAR(irgb)
 110   continue
       do 111 i=143,256
       irgb=float(255*(i-142))/57.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(3,i)=CHAR(irgb)
 111   continue
C-----green
       do 114 i=1,199
       irgb=float(255*(i-29))/56.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(2,i)=CHAR(irgb)
 114   continue
       do 115 i=200,256
       irgb=256.0-float(255*(i-199))/57.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(2,i)=CHAR(irgb)
 115   continue
C------blue
       do 118 i=1,142
       irgb=256.0-float(255*(i-85))/57.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(1,i)=CHAR(irgb)
 118   continue
       do 119 i=143,256
       irgb=0
       cpl(1,i)=CHAR(irgb)
 119   continue
C------last black
       irgb=255
       cpl(1,1)=CHAR(irgb)
       cpl(2,1)=CHAR(irgb)
       cpl(3,1)=CHAR(irgb)
C------first white
       irgb=0
       cpl(1,256)=CHAR(irgb)
       cpl(2,256)=CHAR(irgb)
       cpl(3,256)=CHAR(irgb)
C------extra
       do 120 i=1,256
       irgb=0
       cpl(4,i)=CHAR(irgb)
 120   continue
C--------------------------------end palette---------------
* header 1 (file header ; 1--14 byte)
         headmsw( 1: 2) = 'BM'       ! declaring this is BMP file
         itmp=1024+54+ipix*jpix      ! total file size = header + data
         call num2bit4(itmp,byt4)
         headmsw( 3: 6) = byt4(1:4)
         itmp = 0                    ! may be 0
         call num2bit4(itmp,byt4)
         headmsw( 7: 8) = byt4(1:2)
         headmsw( 9:10) = byt4(1:2)
         itmp = 1078              ! must be 54 : total length of header
         call num2bit4(itmp,byt4)
         headmsw(11:14) = byt4(1:4)
* header 2 (bit-map header ; 13--54 byte)
         itmp = 40                   ! 40 : length of bit-map header
         call num2bit4(itmp,byt4)
         headmsw(15:18) = byt4(1:4)
         itmp = ipix                 ! width
         call num2bit4(itmp,byt4)
         headmsw(19:22) = byt4(1:4)
         itmp = jpix                 ! height
         call num2bit4(itmp,byt4)
         headmsw(23:26) = byt4(1:4)
         itmp = 1                    ! must be 1
         call num2bit4(itmp,byt4)
         headmsw(27:28) = byt4(1:2)
         itmp =  8                   ! 8 : color depth in bit.
         call num2bit4(itmp,byt4)
         headmsw(29:30) = byt4(1:2)
         itmp = 0                    ! 0 : no compression
         call num2bit4(itmp,byt4)
         headmsw(31:34) = byt4(1:4)
         headmsw(35:38) = byt4(1:4)
         headmsw(39:42) = byt4(1:4)
         headmsw(43:46) = byt4(1:4)
         itmp = 256                  ! 256 : num. of colors used
         call num2bit4(itmp,byt4)
         headmsw(47:50) = byt4(1:4)
         itmp = 0                    ! 0 : num. of important color
         call num2bit4(itmp,byt4)
         headmsw(51:54) = byt4(1:4)
* writing header part
         write(mpp,'(a54,$)') headmsw(1:54)
* writing color palette data
         itmp = 256 * 4
         write(frmtstr,'(''('',i8.8,''A,$)'')') itmp
         write(mpp,fmt=frmtstr)
     &      ((cpl(k,i),k=1,4),i=1,256) 
	     return
		 end