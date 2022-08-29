       subroutine num2bit4(inum,byt4)
C***********************************************************************
C --------------------------------------
C convert number to 8-bit characters
C --------------------------------------
       character*4 byt4
       character*1 byta,bytb,bytc,bytd
       itmp1 = inum
       itmp2 = itmp1 / 256**3
       bytd = char(itmp2)
       itmp1 =-itmp2 * 256**3 +itmp1
       itmp2 = itmp1 / 256**2
       bytc = char(itmp2)
       itmp1 =-itmp2 * 256**2 +itmp1
       itmp2 = itmp1 / 256
       bytb = char(itmp2)
       itmp1 =-itmp2 * 256    +itmp1
       byta = char(itmp1)
	   byt4=byta//bytb//bytc//bytd
       return
       end