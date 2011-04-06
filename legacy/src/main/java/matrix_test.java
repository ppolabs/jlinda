public class GlobalMembersMatrix_test
{
	//
	// * Copyright (c) 1999-2009 Delft University of Technology, The Netherlands
	// *
	// * This file is part of Doris, the Delft o-o radar interferometric software.
	// *
	// * Doris program is free software; you can redistribute it and/or modify
	// * it under the terms of the GNU General Public License as published by
	// * the Free Software Foundation; either version 2 of the License, or
	// * (at your option) any later version.
	// *
	// * Doris is distributed in the hope that it will be useful,
	// * but WITHOUT ANY WARRANTY; without even the implied warranty of
	// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	// * GNU General Public License for more details.
	// *
	// * You should have received a copy of the GNU General Public License
	// * along with this program; if not, write to the Free Software
	// * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
	// *
	// * Publications that contain results produced by the Doris software should
	// * contain an acknowledgment. (For example: The interferometric processing
	// * was performed using the freely available Doris software package developed
	// * by the Delft Institute of Earth Observation and Space Systems (DEOS), Delft
	// * University of Technology, or include a reference to: Bert Kampes and
	// * Stefania Usai. \"Doris: The Delft Object-oriented Radar Interferometric
	// * software.\" In: proceedings 2nd ITC ORS symposium, August 1999. (cdrom)).
	// *
	// 
	//***************************************************************
	// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/matrixbk_test.cc,v $      *
	// * $Revision: 3.11 $                                            *
	// * $Date: 2005/10/18 13:46:51 $                                 *
	// * $Author: kampes $                                            *
	// *                                                              *
	// * Test program for matrix class.                               *
	// * Compile with Doris Makefile.                                 *
	// #%// BK 17-Aug-2000                                            *
	// ***************************************************************


	public static bk_messages TRACE;
	public static bk_messages DEBUG;
	public static bk_messages INFO;
	public static bk_messages PROGRESS;
	public static bk_messages WARNING;
	public static bk_messages ERROR;
	public static bk_messages matDEBUG;
	public static bk_messages matERROR;
	// ______ used in matrix_bk.cc to keep track of total allocated ______
	public static uint totalallocated =0; // [B] matrixdebuging





	//***************************************************************
	// *    main                                                      *
	// #%// BK 17-Aug-2000                                            *
	// ***************************************************************
	static int Main(int argc, RefObject<String[]> argv)
	  {
	  final int L = 4; // numlines A
	  final int P = 2; // numcols A
	  int i;
	  int j;

	  matrix<complr4> A = new matrix(L, P);
	  matrix<complr4> B = new matrix(P, L+1);
	  for (i =0; i<A.lines(); ++i)
		for (j =0; j<A.pixels(); ++j)
		  A(i,j) = complr4((float)i/2.3+j,(float)j/2+i);
	  for (i =0; i<B.lines(); ++i)
		for (j =0; j<B.pixels(); ++j)
		  B(i,j) = complr4((float)(j+i)/2.3+i,(float)i/(2.2+i-j *j));

	  // ______ Show matrices for testing ______
	  System.out.print("matrix A(");
	  System.out.printA.lines();
	  System.out.print(",");
	  System.out.printA.pixels();
	  System.out.print("):\n");
	  A.showdata();
	  System.out.print("matrix B(");
	  System.out.printB.lines();
	  System.out.print(",");
	  System.out.printB.pixels();
	  System.out.print("):\n");
	  B.showdata();
	  System.out.print("\n\n---------------------\n\n");

	  // ______ Test - ______
	  B = complr4(2.) * -B;
	  System.out.print("B = 2*-B(");
	  System.out.printB.lines();
	  System.out.print(",");
	  System.out.printB.pixels();
	  System.out.print("):\n");
	  B.showdata();
	  System.out.print("\n\n---------------------\n\n");

	  // ______ Test * ______
	  matrix<complr4> C = A *B;
	  System.out.print("matrix C(");
	  System.out.printC.lines();
	  System.out.print(",");
	  System.out.printC.pixels();
	  System.out.print(") = A*B\n");
	  C.showdata();
	  System.out.print("\n\n---------------------\n\n");

	  // ______ Test matTxmat ______
	  matrix<complr4> C1 = matTxmat(transpose(A), B);
	  System.out.print("matrix C(");
	  System.out.printC1.lines();
	  System.out.print(",");
	  System.out.printC1.pixels();
	  System.out.print(") = A*B (using matTxmat)\n");
	  C1.showdata();
	  System.out.print("\n\n---------------------\n\n");

	  // ______ Test matTxmat ______
	  matrix<complr4> C2 = matxmatT(A, transpose(B));
	  System.out.print("matrix C(");
	  System.out.printC2.lines();
	  System.out.print(",");
	  System.out.printC2.pixels();
	  System.out.print(") = A*B (using matxmatT)\n");
	  C2.showdata();
	  System.out.print("\n\n---------------------\n\n");
	  // --- check result ---
	  System.out.print("checksum1 = ");
	  System.out.printmaxmagnitude(C-C1);
	  System.out.print("\n");
	  System.out.print("checksum2 = ");
	  System.out.printmaxmagnitude(C-C2);
	  System.out.print("\n");



	  // --- Conjugated ---
	  matrix<complr4> Bconj = conj(B);
	  System.out.print("conjugated:\n");
	  Bconj.showdata();
	  System.out.print("\n\n---------------------\n\n");

	  // ______ Test 1d fft ______
	  RefObject<matrix<complr4>> TempRefObject = new RefObject<matrix<complr4>>(A);
	  fft(TempRefObject, 1);
	  A = TempRefObject.argvalue;
	  System.out.print("1d fft over columns of matrix A(");
	  System.out.printA.lines();
	  System.out.print(",");
	  System.out.printA.pixels();
	  System.out.print("):\n");
	  A.showdata();
	  System.out.print("\n\n---------------------\n\n");

	//    // ______ Test operator = ______
	//    for (i=0; i<10; ++i)
	//      {
	//      cout << i << ": A=B\n";
	//      A=B;
	//      cout << i << ": B=B\n";
	//      B=B;
	//      cout << i << ": C=A\n";
	//      matrix<complr4> C=A;
	//      cout << "\n";
	//      }

	  // ______ Test find function ______
	  // matrix<int32> indexNaN = A.find(2);

	  // ______ Test mypow function ______
	  matrix<Float> R = new matrix(L, P);
	  for (i =0; i<R.lines(); ++i)
		for (j =0; j<R.pixels(); ++j)
		  R(i,j) = (float)(i/2.3+j *i *j);

	  System.out.print("\n\n-RRRRRRRRRRRRRRRRRRR-\n\n");
	  R.showdata();
	  R.mypow(1.5);
	  R.showdata();
	// ?? but does B*= A^C work?
	// ?? but does B*= A.mypow(C) work?

	  // ______ Test CR4*R4 function ______
	  System.out.print("\n\n-RRRRRRRRRRRRRRRRRRR-\n\n");
	  matrix<complr4> Q = new matrix(L, P);
	  for (i =0; i<Q.lines(); ++i)
		for (j =0; j<Q.pixels(); ++j)
		  Q(i,j) = complr4((float)i/2.3+j,(float)j/2+i);
	  System.out.print("\n\n- Q.showdata()\n\n");
	  Q.showdata();
	  System.out.print("\n\n- R.showdata()\n\n");
	  R.showdata();
	  Q*=R;
	  System.out.print("\n\n- (Q*=R).showdata()\n\n");
	  Q.showdata();

	  // ______ Test conjugated function ______
	  System.out.print("\n\n- CONJUGATED: Q.conj()\n\n");
	  Q.conj();
	  Q.showdata();
	  System.out.print("\n\n- Q2 = CONJ(Q)\n\n");
	  matrix<complr4> Q2 = conj(Q);
	  Q2.showdata();
	return 0;

	  // ______ Test smooth function ______
	  // R.smooth(1);
	  // ______ Test pow ^ operator ______
	  // R = R^2.0;

	  // ______ Test wshift function ______
	  System.out.print("\n\n- WSHIFT -\n\n");
	  matrix<real4> AA = new matrix(1, 7);
	  for (i =0; i<AA.lines(); ++i)
		for (j =0; j<AA.pixels(); ++j)
		  AA(i,j) = i+j;
	  AA.showdata();
	  System.out.print("wshift AA, -2\n");
	  wshift(AA,-2);
	  AA.showdata();

	  // ______ Test diagxmat R4*CR4 ______
	  matrix<complr4> QQ = new matrix(7, 7);
	  for (i =0; i<QQ.lines(); ++i)
		for (j =0; j<QQ.pixels(); ++j)
		  QQ(i,j) = complr4((float)i/2.3+j,(float)j/2+i);
	  QQ.showdata();
	  matrix<complr4> BB = diagxmat(AA, QQ);
	  System.out.print("BB=diagxmat(AA,QQ)\n");
	  BB.showdata();

	  // ______ Test cos/sin complr4 operator ______
	  System.out.print("VVVVVVVVVVVVVVVVVVVVVVVVVVVVV\n");
	  matrix<real8> r8A = new matrix(3, 5);
	  r8A = 5.3;
	  System.out.print("r8A=5.3:\n");
	  r8A.showdata();

	  for (i =0; i<r8A.lines(); ++i)
		for (j =0; j<r8A.pixels(); ++j)
		  r8A(i,j) = (3.14/180.0)*((0.5+i-j)*100); // rad
	  System.out.print("r8A=radians:\n");
	  r8A.showdata();
	  matrix<real8> cosr8A = Math.cos(r8A);
	  System.out.print("\ncosr8A:\n");
	  cosr8A.showdata();
	  //
	  matrix<real8> r8B = r8A+0.2;
	  matrix<real8> sinr8B = Math.sin(r8B);
	  System.out.print("\nsinr8B:\n");
	  sinr8B.showdata();
	  //
	  matrix<complr4> cr4AB = mat2cr4(r8A,r8B);
	  System.out.print("\ncomplr4(r8A,r8B):\n");
	  cr4AB.showdata();

	  // ______ Test mysort rows on first col ______
	  System.out.print("\n\nmysort(real4) matrix ascending on rows\n");
	  matrix<real4> SS = new matrix(5, 3);
	  SS(0,0)=2.1;
	  SS(0,1)=10.1;
	  SS(0,2)=500.1;
	  SS(1,0)=1.1;
	  SS(1,1)=20.1;
	  SS(1,2)=400.1;
	  SS(2,0)=3.1;
	  SS(2,1)=30.1;
	  SS(2,2)=300.1;
	  SS(3,0)=5.1;
	  SS(3,1)=40.1;
	  SS(3,2)=200.1;
	  SS(4,0)=1.1;
	  SS(4,1)=50.1;
	  SS(4,2)=100.1;
	  System.out.print("original data (SS)\n");
	  SS.showdata();
	  System.out.print("mysort2(SS)\n");
	  RefObject<matrix<real4>> TempRefObject2 = new RefObject<matrix<real4>>(SS);
	  mysort2(TempRefObject2);
	  SS = TempRefObject2.argvalue;
	  SS.showdata();

	  matrix<int32> SSi = new matrix(5, 3);
	  SSi(0,0)=2;
	  SSi(0,1)=10;
	  SSi(0,2)=500;
	  SSi(1,0)=1;
	  SSi(1,1)=20;
	  SSi(1,2)=400;
	  SSi(2,0)=3;
	  SSi(2,1)=30;
	  SSi(2,2)=300;
	  SSi(3,0)=5;
	  SSi(3,1)=40;
	  SSi(3,2)=200;
	  SSi(4,0)=1;
	  SSi(4,1)=50;
	  SSi(4,2)=100;
	  System.out.print("original data (SSi)\n");
	  SSi.showdata();
	  System.out.print("mysort2(SSi)\n");
	  RefObject<matrix<real4>> TempRefObject3 = new RefObject<matrix<real4>>(SSi);
	  mysort2(TempRefObject3);
	  SSi = TempRefObject3.argvalue;
	  SSi.showdata();

	  // ______ Test oversample ______
	  System.out.print("\n\nmysort(real4) matrix ascending on rows\n");
	  matrix<real4> OVS = new matrix(4, 4);
	  OVS(0,0)=2.1;
	  OVS(0,1)=10.1;
	  OVS(0,2)=500.1;
	  OVS(1,0)=1.1;
	  OVS(1,1)=20.1;
	  OVS(1,2)=400.1;
	  OVS(2,0)=3.1;
	  OVS(2,1)=30.1;
	  OVS(2,2)=300.1;
	  OVS(3,0)=5.1;
	  OVS(3,1)=40.1;
	  OVS(3,2)=200.1;
	  System.out.print("original data (OVS)\n");
	  OVS.showdata();
	  System.out.print("mysort2(SS)\n");
	  matrix<real4> OVS_ = oversample(OVS,2,2);
	  OVS_.showdata();

	  // ______ Test fliplr/flipud ______
	  System.out.print("\n\nOVS.showdata()\n");
	  OVS.showdata();
	  System.out.print("OVS.flipud()\n");
	  OVS.flipud();
	  OVS.showdata();
	  System.out.print("OVS.fliplr() (ie both up and lr)\n");
	  OVS.fliplr();
	  OVS.showdata();

	  System.out.print("\n\nNormal termination.\n\n");
	  return 0;
	  } // END main
}


//----------------------------------------------------------------------------------------
//	Copyright © 2006 - 2008 Tangible Software Solutions Inc.
//
//	This class is used to simulate the ability to pass arguments by reference in Java.
//----------------------------------------------------------------------------------------
final class RefObject<T>
{
	T argvalue;
	RefObject(T refarg)
	{
		argvalue = refarg;
	}
}