public class GlobalMembersMatrixbk
{

	// ______ Keep track of total allocated memory ______
	// ______ This will work upto 2^31 bytes (2GB) ______
	// ______ there is a problem here, but how to declare a global? ______
	// ______ if called from different files, not correct bookkeeping ______
	// ______ therefor we define it above main in processor.cc ______
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __DEBUGMAT2
//C++ TO JAVA CONVERTER NOTE: 'extern' variable declarations are not required in Java:
	//  extern uint totalallocated; // [B] extern bk_messages matERROR;
	//#else
//C++ TO JAVA CONVERTER NOTE: 'extern' variable declarations are not required in Java:
	//extern bk_messages matERROR;
	//#endif

	// ______ message objects, global, set in main ______
//C++ TO JAVA CONVERTER NOTE: 'extern' variable declarations are not required in Java:
	//extern bk_messages matDEBUG;
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>




	// ====== Start of class definition ======
	// ====== Private functions ======
	//***************************************************************
	// * allocate                                                     *
	// *    allocate 2d array in major row order, continuously in     *
	// *    memory (required by VECLIB/FFTW). consider rewriting to vector    *
	// *    of vectors of stl, but may not be cont. in mem.           *
	// *    g++ w/o exception handling, ...                           *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::allocate(uint numlines, uint numpixels)
	public void matrix<Type>allocate(uint numlines, uint numpixels) // allocator
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (numlines ==0 || numpixels ==0)
		{
		matERROR << "Allocation impossible: size (l,p): " << numlines << ", " << numpixels;
		matERROR.print();
		}
	//  #endif
	  nrows = numlines;
	  ncols = numpixels;
	  nsize = numlines *numpixels;
	  // Bert Kampes, 07-Apr-2005: try/catch should work by now...
	  try
	      {
		  data = new Type[numlines]; // get memory : make linear array of pointers
	      }
	  catch (bad_alloc)
	      {
		  matERROR << "code 502: first allocation failed, size: " << numlines << ", " << numpixels;
		  matERROR.print();
	      }
	  try
	      {
		  data[0] = new Type[nsize]; // get memory : get the first pointer for the memory block, later we'll fill in all the address
	      }
	  catch(bad_alloc)
	      {
		  matERROR << "code 502: second allocation failed: size: " << numlines << ", " << numpixels;
		  matERROR.print();
	      }
	  for (register uint i =1; i<numlines; i++)
	      data[i] = data[i-1]+numpixels; // start at 0,0

	  //C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	  //  #if __DEBUGMAT2
	  uint allocated = sizeof(Type)*nsize; // [B]
	  totalallocated += allocated; // [B]
	  matDEBUG << "allocated   matrix(" << numlines << "," << numpixels << ") at: " << &data[0][0] << " (" << setw(10) << allocated << "B, total: " << setw(10) << totalallocated << "B)";
	  matDEBUG.print();
	  //  #endif
	  } // END allocate
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>


	//***************************************************************
	// * initialize = allocate + set 0                                *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::initialize(uint numlines, uint numpixels)
	public void matrix<Type>initialize(uint numlines, uint numpixels)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	      matDEBUG.print("initialization matrix.");
	      //  #endif
	      allocate(numlines, numpixels);
	      clean();
	  } // END initialize
    //C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
    //ORIGINAL LINE: template <class Type>



	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __DEBUGMAT1
	//***************************************************************
	// * checkindex (l,p)                                             *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::checkindex(uint line, uint pixel) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public void matrix<Type>checkindex(uint line, uint pixel)
	  {
	  if (int32(line) > int32(nrows)-1 || int32(pixel) > int32(ncols)-1)
		{
		matERROR << "Wrong index (l,p)=(" << line << "," << pixel << "); Matrix(" << nrows << "," << ncols << ") at " << &data[0][0];
		matERROR.print();
		}
	  } // END checkindex template <class Type>
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>
	//#else
	//#endif



	// ====== Public functions ======
	// ====== Constructors ======
	//***************************************************************  \\
	// * matrix<real8> A;                                             *  \\
	// * Bert Kampes, 11-Dec-1998                                     *  \\
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>::matrix()
	public matrix<Type>matrix() // constructor (0 arg)
	  {
	  nrows = 0;
	  ncols = 0;
	  nsize = 0;
	  data = 0; // address of pointer array is set to null
	  } // END constructor
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * matrix<real8> A(3,3);                                        *
	// * Bert Kampes, 11-Dec-1998                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>::matrix(uint lines, uint pixels)
	public matrix<Type>matrix(uint lines, uint pixels)
	  {
	  initialize(lines, pixels); // set to 0
	  } // END constructor
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * matrix<real8> A=B;                                           *
	// * copy constructor; avoids default for dynamical memory:       * 
	// * bitwise copy.                                                *
	// * Bert Kampes, 11-Dec-1998                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>::matrix(const matrix<Type>& A)
	public matrix<Type>matrix(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("copy constructor."); // [MA] ex: copy constructor is called when operator << is used.
	//  #endif
	  if (A.nsize)
		{
		allocate(A.nrows, A.ncols);
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
		memcpy(data[0],A.data[0],nsize *sizeof(Type));
		}
	  else
	  { // [MA] allow for NULL matrix returns when dynamic memory nsize =0
		allocate(A.nrows, A.ncols);
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	#if __DEBUGMAT2
		matDEBUG << "new mtx nsize: " << nsize << " datasize: " << nsize *sizeof(Type) << " dataaddress: " << &data << " vs inputadd: " << A << " addressNULLmemblock: " << data[0];
		matDEBUG.print();
	//	#endif
		//memcpy(data[0],0,nsize*sizeof(Type)); // clean() // when nsize=0 and A.data[0] doesn't exist (data=0), such as a null mtx, thus memcpy crashes with "Caught SIGSEGV: Segmentation fault."
											  // this is necessary otherwise it points to a memory block
	   data[0]=0;
	   }
	  } // END constructor
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * matrix<int32> A(win,B)                                       *
	// * Constructor (part)                                           *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// *    Bert Kampes, 31-Mar-1999 memcpy ipv. for_j                *
	// * window starts at 0 (matrix index)                            *
	// #%// BK 25-Oct-2000                                            *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>::matrix (window win, const matrix<Type>& A)
	public matrix<Type>matrix(window win, matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("constructor as part.");
	//  #endif
	  // ______ Check arguments ______
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (win.linehi<win.linelo)
		matERROR.print("constructor (4uint,matrix): win.linehi<linelo ?");
	  if (win.pixhi<win.pixlo)
		matERROR.print("constructor (4uint,matrix): win.pixhi<pixlo ?");
	  A.checkindex(win.linehi, win.pixhi);
	//  #endif
	  // ______ Allocate new matrix and fill ______
	  final uint numlin = win.lines();
	  final uint numpix = win.pixels();
	  final uint sizelin = numpix *sizeof(Type);
	  allocate(numlin, numpix);
	  for(register uint i =0; i<numlin; i++)
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
		memcpy(data[i],A[win.linelo+i]+win.pixlo,sizelin);
	  } // END constructor
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	// ======Destructor======
	//***************************************************************
	// * Destructor                                                   *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>::~matrix()
	public matrix<Type>~matrix()
	  {
	  if (data==0)
		  return;
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  uint deallocated = sizeof(Type)*nsize; // [B]
	  totalallocated -= deallocated; // [B]
	  matDEBUG << "deallocated matrix(" << nrows << "," << ncols << ") at: " << data[0] << " (" << setw(10) << deallocated << "B; total: " << setw(10) << totalallocated << "B)";
	  matDEBUG.print();

	//  #endif
	  data[0] = null; // deallocate
	  data = null; // deallocate
	  data=0; // set to null pointer
	  } // END destructor
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	// ======Data functions======
	//***************************************************************
	// * A.setdata(w)                                                 *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::setdata(Type w)
	public void matrix<Type>setdata(Type w) // sets matrix to constant
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("A.setdata(w)");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pnt = data[0];
	  for (register uint i =0;i<nsize;i++)
		(*pnt++) = Type(w);
	  } // END setdata
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * setdata(i,j,A)                                               *
	// *    put matrix A on l1,p1                                     *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::setdata(uint l1, uint p1, const matrix<Type>& A)
	public void matrix<Type>setdata(uint l1, uint p1, matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  checkindex(l1+A.nrows-1, p1+A.ncols-1); // check far most corner const uint sizelin = A.ncols *sizeof(Type);
	//  #else
	  final uint sizelin = A.ncols *sizeof(Type);
	//  #endif
	  for (register uint i =0;i<A.nrows;i++)
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
		memcpy(data[i+l1]+p1,A[i],sizelin);
	  } // END setdata
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * B.setdata(win, A, winA):                                     *
	// *  set win of B to winA of A                                   *
	// * if win==0 defaults to totalB, winA==0 defaults to totalA     *
	// * first line matrix =0 (?)
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::setdata(window winin, const matrix<Type> &A, window winA)
	public void matrix<Type>setdata(window winin, matrix<Type> A, window winA)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("setdata (win,A,win)");

	//  #endif
	  // ______Check default request______
	  if (winin.linehi == 0 && winin.pixhi == 0)
		{
			winin.linehi = nrows -1;
		 winin.pixhi = ncols -1;
		}
	  if (winA.linehi == 0 && winA.pixhi == 0)
		{
			winA.linehi = A.lines() -1;
		 winA.pixhi = A.pixels() -1;
		}
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __DEBUGMAT1
	  if (((winin.linehi - winin.linelo) != (winA.linehi - winA.linelo)) || ((winin.pixhi - winin.pixlo) != (winA.pixhi - winA.pixlo)))
		matERROR.print("code 901: wrong input.");
	  if (winin.linehi<winin.linelo || winin.pixhi<winin.pixlo)
		matERROR.print("code 901: wrong input.1");
	  if ((winin.linehi > nrows -1) || (winin.pixhi > ncols -1))
		matERROR.print("code 901: wrong input.2");
	  if ((winA.linehi > A.lines() -1) || (winA.pixhi > A.pixels() -1))
		matERROR.print("code 901: wrong input.3");
	//#endif
	  // ______ Fill data ______
	  final uint sizelin = (winA.pixhi - winA.pixlo + 1)*sizeof(Type);
	  for(register uint i =winin.linelo; i<=winin.linehi;i++)
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
		memcpy(data[i]+winin.pixlo,A[i-winin.linelo+winA.linelo]+winA.pixlo,sizelin);
	  } // END setdata
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * B.setdata(A, winA):                                          *
	// *    set total of B to winA of A                               *
	// *    Bert Kampes, 17-Mar-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::setdata(const matrix<Type> &A, window winA)
	public void matrix<Type>setdata(matrix<Type> A, window winA)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("setdata (A,win)");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (((nrows -1) != (winA.linehi - winA.linelo)) || ((ncols -1) != (winA.pixhi - winA.pixlo)))
		matERROR.print("code 901: wrong input.");
	  if ((winA.linehi > A.lines() -1) || (winA.pixhi > A.pixels() -1))
		matERROR.print("code 901: wrong input.3");
	//  #endif
	  // ______ Fill data ______
	  final uint sizelin = ncols *sizeof(Type);
	  for(register uint i =0; i<nrows;i++)
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
		memcpy(data[i],A[i+winA.linelo]+winA.pixlo,sizelin);
	  } // END setdata
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * A=B.getdata(win)                                             *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type> matrix<Type>::getdata(window win) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public matrix<Type> matrix<Type>getdata(window win)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("getdata.");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (win.linehi<win.linelo || win.pixhi<win.pixlo)
		matERROR.print("code 501: matrix::getdata (win): arguments are wrong, l1<l2,p1<p2");
	  checkindex(win.linehi, win.pixhi);
	//  #endif
	  final uint numlin = win.lines();
	  final uint numpix = win.pixels();
	  matrix<Type> Result = new matrix(numlin, numpix); // =data(;
	  for(register uint i =0; i<numlin; i++)
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
		memcpy(Result[i],data[i+win.linelo]+win.pixlo,numpix *sizeof(Type));
	  return Result;
	  } // END getdata
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * A=B.getrow(row)                                              *
	// * rows: 0 1 2 ..                                               *       
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type> matrix<Type>::getrow(uint line) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public matrix<Type> matrix<Type>getrow(uint line)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  checkindex(line, 0);
	//  #endif
	  matrix<Type> Result = new matrix(1, ncols);
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
	  memcpy(Result[0],data[line],ncols *sizeof(Type));
	  return Result;
	  } // END getrow
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * A=B.getcolumn(col)                                           *
	// * cols: 0 1 2 ..                                               *       
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type> matrix<Type>::getcolumn(uint pixel) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public matrix<Type> matrix<Type>getcolumn(uint pixel)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
		checkindex(0, pixel);
	//  #endif
	  matrix<Type> Result = new matrix(nrows, 1);
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntA = data[0]+pixel;
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntR = Result[0];
	  for (register uint i =0; i<nrows; ++i)
		{
		(*pntR++) = *pntA;
		pntA += ncols;
		}
	  return Result;
	  } // END getcolumn
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * B.showdata()                                                 *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// *    Mahmut Arikan, 22-May-2009 - Adjustment to see more digits*
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::showdata() const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public void matrix<Type>showdata() // show all data in matrix
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
		matDEBUG.print("showdata.");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  if (nrows>100 || ncols>15)
		{
		matDEBUG << "matrix (" << nrows << "," << ncols << "); only showing data (0:99,0:9).";
		matDEBUG.print();
		}
	//  #endif
	  final uint L = (nrows<=100) ? nrows : 100;
	  final uint P = (ncols<=10) ? ncols : 10;
	  matDEBUG.precision(11); // [MA] 9 --> 11
	  matDEBUG.width(12); // 10 --> 12
	  //matDEBUG.setf(ios::left, ios::adjustfield); // [MA] this failed by 4.01 version of messages.hh
	  for (register uint i =0; i<L; i++)
		{
		for (register uint j =0; j<P; j++)
		  {
		  // matDEBUG << data[i][j] << " ";
		  matDEBUG << left << data[i][j] << " ";
		  }
		matDEBUG.print(); // prevent too much in ostream
		}
	  matDEBUG.reset();
	  } // END showdata
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * bool v = A.isvector()                                        *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: inline boolean matrix<Type>::isvector() const
//C++ TO JAVA CONVERTER NOTE: Inline methods are not available in Java:
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public boolean matrix<Type>isvector()
	  {
	  return (nrows==1 || ncols==1) ? true : false;
	  } // END isvector
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * uint l = A.lines()                                           *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: inline uint matrix<Type>::lines() const
//C++ TO JAVA CONVERTER NOTE: Inline methods are not available in Java:
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public uint matrix<Type>lines() // return number of lines
	  {
	  return nrows;
	  } // END lines
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * uint p = A.pixels()                                          *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: inline uint matrix<Type>::pixels() const
//C++ TO JAVA CONVERTER NOTE: Inline methods are not available in Java:
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public uint matrix<Type>pixels() // return number of pixels
	  {
	  return ncols;
	  } // END pixels
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * uint s = A.size()                                            *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: inline uint matrix<Type>::size() const
//C++ TO JAVA CONVERTER NOTE: Inline methods are not available in Java:
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public uint matrix<Type>size() // return nsize
	  {
	  return nsize;
	  } // END size
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * A.resize(l,p)                                                *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::resize(uint l1, uint p1)
	public void matrix<Type>resize(uint l1, uint p1)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("resize.");
	//  #endif
	  if (l1 == nrows && p1 == ncols)
		  return;
	  else if (data!=0) // check for allocated memory
		{
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	#if __DEBUGMAT2
		uint deallocated = sizeof(Type)*nsize; // [B]
		totalallocated -= deallocated; // [B]
		matDEBUG << "deallocated matrix(" << nrows << "," << ncols << ") at: " << data[0] << " (" << setw(10) << deallocated << "B; total: " << setw(10) << totalallocated << "B)";
		matDEBUG.print();
	//	#endif
		data[0] = null;
		data = null;
		data=0;
		}
	  initialize(l1, p1); // set to 0
	  } // END resize
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * A.clean()                                                    *
	// *    Bert Kampes, 01-Feb-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::clean()
	public void matrix<Type>clean() // sets 2 zero
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("mtx clean.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memset' has no equivalent in Java:
	  memset(data[0],0,nsize *sizeof(Type));
	  } // END clean
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * B.setrow(row, data)                                          *
	// * rows: 0 1 2 .., should fit exactly                           *       
	// * orientation of vector is disregarded.                        *
	// *    Bert Kampes, 12-Oct-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::setrow(uint line, const matrix<Type> &LINE)
	public void matrix<Type>setrow(uint line, matrix<Type> LINE)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  checkindex(line, 0);
	  if (!(LINE.nrows == 1 || LINE.ncols == 1))
		matERROR.print("setrow: only vector input.");
	  if (LINE.nsize != ncols)
		matERROR.print("setrow: sizeofvector should be same as matrix.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
	  memcpy(data[line],LINE[0],ncols *sizeof(Type));
	  } // END setrow
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * B.setrow(row, scalar)                                        *
	// * rows: 0 1 2 .., should fit exactly                           *       
	// * orientation of vector is disregarded.                        *
	// *    Bert Kampes, 12-Oct-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::setrow(uint line, Type scalar)
	public void matrix<Type>setrow(uint line, Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  checkindex(line, 0);
	//  #endif
	  for (register Type *pntB = &data[line][0]; pntB<=&data[line][ncols-1]; pntB++)
		*pntB = scalar;
	  } // END setrow
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * B.setcolumn(col, COL)                                        *
	// * cols: 0 1 2 ..                                               *       
	// * orientation of vector is disregarded.                        *
	// *    Bert Kampes, 12-Oct-1999                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::setcolumn(uint pixel, const matrix<Type> &COLUMN)
	public void matrix<Type>setcolumn(uint pixel, matrix<Type> COLUMN)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  checkindex(0, pixel);
	  if (!(COLUMN.nrows == 1 || COLUMN.ncols == 1))
		matERROR.print("setcolumn: only vector input.");
	  if (COLUMN.nsize != nrows)
		matERROR.print("setcolumn: sizeofvector should be same as matrix.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntCOL = COLUMN[0];
	  for (register Type *pntB = &data[0][pixel]; pntB<=&data[nrows-1][pixel]; pntB+=ncols)
		*pntB = *pntCOL++;
	  } // END setcolumn
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * B.setcolumn(col, scalar)                                     *
	// * cols: 0 1 2 ..                                               *       
	// #%// BK 25-Sep-2000                                            *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::setcolumn(uint pixel, Type scalar)
	public void matrix<Type>setcolumn(uint pixel, Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  checkindex(0, pixel);
	//  #endif
	  for (register Type *pntB = &data[0][pixel]; pntB<=&data[nrows-1][pixel]; pntB+=ncols)
		*pntB = scalar;
	  } // END setcolumn
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * B.fliplr()                                                   *
	// * Mirror in center vertical (flip left right).                 *
	// *    Bert Kampes, 23-Mar-2000                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::fliplr()
	public void matrix<Type>fliplr()
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("fliplr.");
	//  #endif
	  if (nrows==1)
		{
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		Type *pnt1 = data[0]; // first one
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		Type *pnt2 = data[0]+ncols-1; // last one
		Type tmp = *pnt1;
		for (register int32 i =0; i<int32(ncols/2); ++i) // floor
		  {
		  (*pnt1++) = *pnt2;
		  (*pnt2--) = tmp;
		  tmp = *pnt1;
		  }
		}
	  else
		{
		for (register int32 i =0; i<int32(ncols/2); ++i) // floor
		  {
		  matrix<Type> tmp1 = getcolumn(i);
		  matrix<Type> tmp2 = getcolumn(ncols-i-1);
		  setcolumn(i, tmp2);
		  setcolumn(ncols-i-1, tmp1);
		  }
		}
	  } // END fliplr
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * B.flipud()                                                   *
	// * Mirror in center vertical (flip left right).                 *
	// *    Actually move data around, not pointers, to be sure       *
	// *    veclib works ok, data is cont. in memory.                 *
	// *    Bert Kampes, 23-Mar-2000                                  *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::flipud()
	public void matrix<Type>flipud()
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("flipud.");
	//  #endif
	  if (ncols==1)
		{
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		Type *pnt1 = data[0]; // first one
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		Type *pnt2 = data[0]+nrows-1; // last one
		Type tmp = *pnt1;
		for (register int32 i =0; i<int32(nrows/2); ++i) // floor
		  {
		  (*pnt1++) = *pnt2;
		  (*pnt2--) = tmp;
		  tmp = *pnt1;
		  }
		}
	  else
		{
		for (register int32 i =0; i<int32(ncols/2); ++i) // floor
		  {
		  matrix<Type> tmp1 = getrow(i);
		  matrix<Type> tmp2 = getrow(nrows-i-1);
		  setrow(i, tmp2);
		  setrow(nrows-i-1, tmp1);
		  }
		}
	  } // END flipud
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>







	//***************************************************************
	// ****************************************************************
	// * OVERLOADED OPS                                               *
	// ****************************************************************
	// ***************************************************************




	//***************************************************************
	// * a = A[5][2];                                                 *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: inline Type* matrix<Type>::operator [] (uint line) const
//C++ TO JAVA CONVERTER NOTE: Inline methods are not available in Java:
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
	public <Type> Type matrix<Type>getDefault(uint line)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  checkindex(line, 0);
	//  #endif
	  return data[line];
	  } // END []
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * a = A(i,j);                                                  *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: inline Type& matrix<Type>::operator ()(uint line, uint pixel) const
//C++ TO JAVA CONVERTER NOTE: Inline methods are not available in Java:
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	Type matrix<Type>operator ()(uint line, uint pixel)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  checkindex(line, pixel);
	//  #endif
	  return data[line][pixel];
	  } // END ()
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * matrix<T> B = A(window);                                     *
	// * may not be to efficient cause twice allocation?              *
	// * Bert Kampes, 31-Mar-2000                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type> matrix<Type>::operator ()(const window &win) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator ()(window win)
	  {
	  matrix<Type> Res = new matrix(win, this);
	  return Res;
	  } // END (win)
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * matrix<T> B = A(uint,uint,uint,uint);                        *
	// * Bert Kampes, 31-Mar-2000                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type> matrix<Type>::operator ()(const uint &l0, const uint &lN, const uint &p0, const uint &pN) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator ()(uint l0, uint lN, uint p0, uint pN)
	  {
	  final window win = new window(l0,lN,p0,pN);
	  matrix<Type> Res = new matrix(win, this);
	  return Res;
	  } // END (4 uint)
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// *  =                                                           *
	// * Bert Kampes, 14-Jan-1999                                     *
	// * Mahmut Arikan, 19-May-2009 Null matrix  assignment           *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator = (const matrix<Type>& A)
//C++ TO JAVA CONVERTER NOTE: This 'CopyFrom' method was converted from the original C++ copy assignment operator:
	public matrix<Type> matrix<Type>CopyFrom(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("operator =");
	//  #endif
	  //cerr << "operator =\n";
	  if (this != A) // prevent copy to itself
		{
		if (A.nsize) // if allocated : disable this to return null matrix
		  {
		  if (data != 0) // if allocated
			{
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//		#if __DEBUGMAT2
			uint deallocated = sizeof(Type)*nsize; // [B]
			totalallocated -= deallocated; // [B]
			matDEBUG << "deallocated matrix(" << nrows << "," << ncols << ") at: " << data[0] << " (" << setw(10) << deallocated << "B; total: " << setw(10) << totalallocated << "B)";
			matDEBUG.print();
	//		#endif
			data[0] = null;
			data = null;
			data = 0;
			}
		  allocate(A.nrows, A.ncols);
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
		  memcpy(data[0],A.data[0],nsize *sizeof(Type));
		  }
		else
		  { // [MA] allow for NULL matrix returns when dynamic memory nsize =0
			cerr << "op =    : NULL matrix is assigned.";
			allocate(A.nrows, A.ncols);
			data[0]=0; // point to no memory block
		   }
		}
	  return this;
	  } // END =
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>


	//***************************************************************
	// *  =                                                           *
	// #%// BK 09-Nov-2000                                            *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator = (const Type scalar)
//C++ TO JAVA CONVERTER NOTE: This 'CopyFrom' method was converted from the original C++ copy assignment operator:
	public matrix<Type> matrix<Type>CopyFrom(Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("operator = (scalar)");
	//  #endif
	  setdata(scalar);
	  return this;
	  } // END = (scalar)
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>




	//***************************************************************
	// * C *= 5.0;                                                    *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator *= (Type scalar)
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator *=(Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("*=");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (!nsize)
		matERROR.print("matrix:: *= with empty matrix.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntmat = data[0];
	  for (register uint i =0; i<nsize; ++i)
		(*pntmat++) *= scalar;
	  return this;
	  } // END *= scalar
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>




	//***************************************************************
	// * C *= A;      pointwise multiplication, a,c same size         *
	// * Bert Kampes, 06-Oct-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator *= (const matrix<Type> &A)
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator *=(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("*= pointwise");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (nrows != A.lines() || ncols != A.pixels())
		matERROR.print("matrix:: *= matrices must be same size.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntmat = data[0];
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntA = A[0];
	  //for (register uint i=0; i<nsize; i++)
	  //  (*pntmat++) *= (*pntA++); 
	  // changed by FvL (for g++/gcc > 4.0):
	  for (register uint i =0; i<nsize; i++)
		{
		(*pntmat++) *= (*pntA);
		*pntA++;
		}
	  return this;
	  } // END *= matrices
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>




	//***************************************************************
	// * C /= 5.0;                                                    *
	// * Bert Kampes, 12-Oct-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator /= (Type scalar)
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator /=(Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("/= scalar");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (!nsize)
		matERROR.print("matrix:: /= with empty matrix.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		Type *pntmat = data[0];
		for (register uint i =0; i<nsize; i++)
		  (*pntmat++) /= scalar;
		return this;
	  } // END /= scalar
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>




	//***************************************************************
	// * C /= A;      pointwise division, a,c same size               *
	// * Bert Kampes, 06-Oct-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator /= (const matrix<Type> &A)
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator /=(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("/=");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (nrows != A.lines() || ncols != A.pixels())
		matERROR.print("matrix:: /= matrices must be same size.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntmat = data[0];
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntA = A[0];
	  //for (register uint i=0; i<nsize; i++)
	  //  (*pntmat++) /= (*pntA++); 
	  // changed by FvL (for g++/gcc > 4.0):
	  for (register uint i =0; i<nsize; i++)
		{
		(*pntmat++) /= (*pntA);
		*pntA++;
		}
	  return this;
	  } // END /= matrices
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C -= A;                                                      *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator -= (const matrix<Type>& A)
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator -=(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("-= mat");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (nrows != A.nrows || ncols != A.ncols)
		matERROR.print("error dimensions.");
	  if (!A.nsize)
		matERROR.print("matrix:: -= with empty matrices.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntmat = data[0];
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntA = A.data[0];
	  //for (register uint i=0; i<nsize; i++)
	  //  (*pntmat++) -= (*pntA++); 
	  // changed by FvL (for g++/gcc > 4.0):
	  for (register uint i =0; i<nsize; i++)
		{
		(*pntmat++) -= (*pntA);
		*pntA++;
		}
	  return this;
	  } // END -=
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C -= 5.0;                                                    *
	// * Bert Kampes, 26-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator -= (Type scalar)
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator -=(Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("-= scalar");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (!nsize)
		matERROR.print("matrix:: -= with empty matrix.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntmat = data[0];
	  for (register uint i =0; i<nsize; i++)
		(*pntmat++) -= scalar;
	  return this;
	  } // END -= scalar
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>




	//***************************************************************
	// * C += A;                                                      *
	// * Bert Kampes, 14-Jan-1999                                     *
	// #%// Bert Kampes, 10-Apr-2005  (why const?)
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator += (const matrix<Type>& A)
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator +=(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("+=");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (nrows != A.nrows || ncols != A.ncols)
		matERROR.print("error dimensions.");
	  if (!A.nsize)
		matERROR.print("matrix:: += with empty matrices.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntmat = data[0];
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntA = A.data[0];
	  //for (register uint i=0; i<nsize; i++)
	  //  (*pntmat++) += (*pntA++); 
	  // changed by FvL (for g++/gcc > 4.0):
	  for (register uint i =0; i<nsize; i++)
		{
		(*pntmat++) += (*pntA);
		*pntA++;
		}
	  return this;
	  } // END +=
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>




	//***************************************************************
	// * C += 5.0;                                                    *
	// * Bert Kampes, 26-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: matrix<Type>& matrix<Type>::operator += (Type scalar)
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	matrix<Type> matrix<Type>operator +=(Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("+= scalar");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (!nsize)
		matERROR.print("matrix:: += with empty matrix.");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntmat = data[0];
	  for (register uint i =0; i<nsize; i++)
		(*pntmat++) += scalar;
	  return this;
	  } // END += scalar
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * A.conj();complex conjugated                                  *
	// * Bert Kampes, 18-Oct-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::conj()
	public void matrix<Type>conj()
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("conj");
	//  #endif
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntmat = data[0];
	  for (register uint i =0; i<nsize; ++i)
		{
		(*pntmat) = Type(pntmat.real(), -pntmat.imag());
		pntmat++;
		}
	  } // END conj()
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * if (A==scalar) all elements equal scalar                     *
	// * Bert Kampes, 08-Oct-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: boolean matrix<Type>::operator == (Type scalar) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	boolean matrix<Type>operator ==(Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("== (scalar)");
	//  #endif
	  if (nsize == 0)
		return false;
	  boolean same = true;
	  Type pnt = data[0];
	  for (register int32 i =0; i<nsize; ++i)
		{
		if (( pnt) != scalar)
		  {
		  same = false;
		  break;
		  }
		}
	  return same;
	  } // END ==
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * if (A==B)                                                    *
	// * Bert Kampes, 08-Oct-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: boolean matrix<Type>::operator == (const matrix<Type> &A) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	boolean matrix<Type>operator ==(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __DEBUGMAT2
	  matDEBUG.print("==");
	//#endif
	  if ((A.lines() != nrows) || (A.pixels() != ncols))
		return false;
	  boolean same = true;
	  Type pnt = data[0];
	  Type pntA = A[0];
	  for (register uint i =0; i<nsize; ++i)
		{
		if (( pnt) != ( pntA))
		  {
		  same = false;
		  break;
		  }
		}
	  return same;
	  } // END ==
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * if (A!=scalar)                                               *
	// * Bert Kampes, 08-Oct-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: boolean matrix<Type>::operator != (Type scalar) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	boolean matrix<Type>operator !=(Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __DEBUGMAT2
	  matDEBUG.print("!= (scalar)");
	//#endif
	  if ( this == scalar)
		return false;
	  else
		return true;
	  } // END !=
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * if (A!=B)                                                    *
	// * Bert Kampes, 08-Oct-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: boolean matrix<Type>::operator != (const matrix<Type> &A) const
//C++ TO JAVA CONVERTER WARNING: 'const' methods are not available in Java:
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	boolean matrix<Type>operator !=(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//#if __DEBUGMAT2
	  matDEBUG.print("!=");
	//#endif
	  if ( this == A)
		return false;
	  else
		return true;
	  } // END !=
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>




	// ++++++++++++++++++++++++++
	// template functions, stupid place, but else not found?
	// specializations in matrixspecs.cc are NOT used before these?
	// BK 16-Aug-2000
	// ++++++++++++++++++++++++++
	//***************************************************************
	// * C = A * B;                                                   *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator *(matrix<Type> A, matrix<Type> B)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("matrices * (no veclib)");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (A.pixels() != B.lines())
		matERROR.print("matrix::operator *: multiplication not possible");
	  if (!A.size())
		matERROR.print("matrix:: operator * with empty matrices.");
	//  #endif
	  // ______ Straightforward, no veclib _____
	  // ______ it is probably worth to make this faster, e.g.,
	  // ______ using pointers (using A[i][j] is same speed)
	  // ______ and not initializing Result
	  matrix<Type> Result = new matrix(A.lines(), B.pixels());
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register Type sum = Type(0.0);
	  Type sum = new Type(0.0);
	  // --- use straightforward notation for slow matrix access --------
	//  #define NO_POINTERS
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if NO_POINTERS
	  for (register uint i =0; i<Result.lines(); i++)
		{
		for (register uint j =0; j<Result.pixels(); j++)
		  {
		  for (register uint k =0; k<A.pixels(); k++)
			{
			sum += A(i,k) * B(k,j);
			}
		  Result(i,j) = sum;
		  sum = Type(0.0); // complex requires this
		  }
		}
	  // --- use pointers for faster matrix access --------------------
	  // Bert Kampes, 13-Oct-2005: tested OK.
	//  #else
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntR = Result[0]; // use to point from first to last element
	  for (register uint i =0; i<Result.lines(); i++)
		{
		for (register uint j =0; j<Result.pixels(); j++)
		  {
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		  Type *pntA = A[i]; // point to first element this row of A
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		  Type *pntB = B[0]+j; // point to first element in column of B
		  for (register uint k =0; k<A.pixels(); k++)
			{
			//sum  += (*pntA++) * (*pntB); // pointer over row A
			// changed by FvL (for g++/gcc > 4.0):
			sum += (*pntA) * (*pntB); // pointer over row A
			*pntA++;
			pntB += B.pixels(); // point to next element in this column of B
			}
		  (*pntR++) = sum; // matrix is linear in memory, major-row
		  sum = Type(0.0); // complex requires this
		  }
		}
	//  #endif
	  return Result;
	  } // END *
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = 5.0 * B;                                                 *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator *(matrix<Type> A, Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("scalar *");
	//  #endif
	  // ______Perform multiplication______
	  matrix<Type> Result =A;
	  return Result *= scalar; // checks in *=
	  } // END * scalar
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = B * 5.0;                                                 *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator *(Type scalar, matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("* scalar.");
	//  #endif
	  return A *scalar; // checks in *=
	  } // END scalar *
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = A / 5;                                                   *
	// * Bert Kampes, 04-Apr-2000                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator /(matrix<Type> A, Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("/");
	//  #endif
	  matrix<Type> Result = A;
	  return Result /= scalar; // checks are performed here.
	  } // END /
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = A / B;                                                   *
	// * Bert Kampes, 04-Apr-2000                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator /(matrix<Type> A, matrix<Type> B)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("/");
	//  #endif
	  matrix<Type> Result = A;
	  return Result /= B; // checks are performed here.
	  } // END /
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = A - B;                                                   *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator -(matrix<Type> A, matrix<Type> B)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("-");
	//  #endif
	  matrix<Type> Result = A;
	  return Result -= B; // checks are performed here.
	  } // END - (binary)
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = A - 5;                                                   *
	// * Bert Kampes, 04-Apr-2000                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator -(matrix<Type> A, Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("-");
	//  #endif
	  matrix<Type> Result = A;
	  return Result -= scalar; // checks are performed here.
	  } // END - (binary)
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = A + B;                                                   *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator +(matrix<Type> A, matrix<Type> B)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("+");
	//  #endif
	  matrix<Type> Result = B;
	  return Result += A; // checks are in +=
	  } // END + (binary)
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = A + 5;                                                   *
	// * Bert Kampes, 04-Apr-2000                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator +(matrix<Type> A, Type scalar)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("+");
	//  #endif
	  matrix<Type> Result = A;
	  return Result += scalar; // checks are performed here.
	  } // END + (binary)
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * a = max(A)                                                   *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
	public static <Type> Type max(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("max");
	//  #endif
	  Type m =A(0,0);
	  for (register uint i =0; i<A.lines(); ++i)
		for (register uint j =0; j<A.pixels(); ++j)
		  if (A(i,j)>m)
			  m =A(i,j);
	  return m;
	  } // END max
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * a = max(A,linemax,pixelmax)                                  *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
	public static <Type> Type max(matrix<Type> A, RefObject<uint> line, RefObject<uint> pixel)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("max");
	//  #endif
	  Type m =A(0,0);
	  for (register uint i =0; i<A.lines(); ++i)
		for (register uint j =0; j<A.pixels(); ++j)
		  if (A(i,j)>=m)
			{
			m = A(i,j);
			line.argvalue = i;
			pixel.argvalue = j;
			}
	  return m;
	  } // END max
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * a = min(A)                                                   *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
	public static <Type> Type min(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("min");
	//  #endif
	  Type m =A(0,0);
	  for (register uint i =0; i<A.lines(); ++i)
		for (register uint j =0; j<A.pixels(); ++j)
		  if (A(i,j)<m)
			  m =A(i,j);
	  return m;
	  } // END min
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * a = min(A,linemax,pixelmax)                                  *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
	public static <Type> Type min(matrix<Type> A, RefObject<uint> line, RefObject<uint> pixel)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("min");
	//  #endif
	  Type m =A(0,0);
	  for (register int32 i =0; i<A.lines(); i++)
		for (register int32 j =0; j<A.pixels(); j++)
		  if (A(i,j)<=m)
			{
			m = A(i,j);
			line.argvalue = i;
			pixel.argvalue = j;
			}
	  return m;
	  } // END min
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C=matTxmat(A,B) C=trans(A)*B; specialized for veclib         *
	// *    Bert Kampes, 22-Feb-1999                                  *
	// ***************************************************************
	public static <Type> matrix<Type> matTxmat(matrix<Type> A, matrix<Type> B)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
		matDEBUG.print("matTxmat: no veclib");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (A.lines() != B.lines())
		matERROR.print("matTxmat: size A,B: input is A,B; computed is trans(A)*B.");
	//  #endif
	  matrix<Type> Result = new matrix(A.pixels(), B.pixels());
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register Type sum = Type(0.0);
	  Type sum = new Type(0.0);
	  // --- use straightforward notation for slow matrix access --------
	//  #define NO_POINTERS
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if NO_POINTERS
	  for (register uint i =0; i<Result.lines(); i++)
		{
		for (register uint j =0; j<Result.pixels(); j++)
		  {
		  for (register uint k =0; k<A.lines(); k++)
			{
			sum += A(k,i) * B(k,j);
			}
		  Result(i,j) = sum;
		  sum = Type(0.0);
		  }
		}
	  // --- use pointers for faster matrix access --------------------
	//  #else
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntR = Result[0]; // use to point from first to last element
	  for (register uint i =0; i<Result.lines(); i++)
		{
		for (register uint j =0; j<Result.pixels(); j++)
		  {
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		  Type *pntA = A[0]+i; // point to first element this row of A
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		  Type *pntB = B[0]+j; // point to first element in column of B
		  for (register uint k =0; k<A.lines(); k++)
			{
			sum += (*pntA) * (*pntB); // pointer over row A
			pntA += A.pixels();
			pntB += B.pixels(); // point to next element in this column of B
			}
		  (*pntR++) = sum; // matrix is linear in memory, major-row
		  sum = Type(0.0); // complex requires this
		  }
		}
	//  #endif
	  return Result;
	  } // END matTxmat
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C=matxmatT(A,B) C=A*trans(B); specialized for veclib         *
	// *    Bert Kampes, 22-Feb-1999                                  *
	// ***************************************************************
	public static <Type> matrix<Type> matxmatT(matrix<Type> A, matrix<Type> B)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
		matDEBUG.print("matxmatT: no veclib");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (A.pixels() != B.pixels())
		matERROR.print("matxmatT: size A,B: input is A,B; computed is A*trans(B).");
	//  #endif
//C++ TO JAVA CONVERTER NOTE: 'register' variable declarations are not supported in Java:
//ORIGINAL LINE: register Type sum = Type(0.0);
	  Type sum = new Type(0.0);
	  matrix<Type> Result = new matrix(A.lines(), B.lines());
	  // --- use straightforward notation for slow matrix access --------
	//  #define NO_POINTERS
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if NO_POINTERS
	  for (register uint i =0; i<Result.lines(); i++)
		{
		for (register uint j =0; j<Result.pixels(); j++)
		  {
		  for (register uint k =0; k<A.pixels(); k++)
			{
			sum += A(i,k) * B(j,k);
			}
		  Result(i,j) = sum;
		  sum = Type(0.0);
		  }
		}
	  // --- use pointers for faster matrix access --------------------
	//  #else
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntR = Result[0]; // use to point from first to last element
	  for (register uint i =0; i<Result.lines(); i++)
		{
		for (register uint j =0; j<Result.pixels(); j++)
		  {
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		  Type *pntA = A[i]; // point to first element this row of A
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		  Type *pntB = B[j]; // point to first element in column of B
		  for (register uint k =0; k<A.pixels(); k++)
			{
			//sum  += (*pntA++) * (*pntB++); // pointer over row A
			// changed by FvL (for g++/gcc > 4.0):
			sum += (*pntA) * (*pntB); // pointer over row A
			*pntA++;
			*pntB++;
			}
		  (*pntR++) = sum; // matrix is linear in memory, major-row
		  sum = Type(0.0); // complex requires this
		  }
		}
	//  #endif
	  return Result;
	  } // END matxmatT
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * dumpasc(file,A);                                             *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
	public static <Type> void dumpasc(String file, matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("dumpasc to file");
	//  #endif
	  ofstream fo = new ofstream(file,ios.out | ios.trunc);
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __LINE__ macro:
	//C++ TO JAVA CONVERTER TODO TASK: There is no direct equivalent in Java to the C++ __FILE__ macro:
	  matassert(fo,file,__FILE__,__LINE__);
	  fo.precision(3);
	  fo.width(11);
	  fo.setf(ios.fixed);
	  for (register int32 i =0; i<A.lines(); ++i)
		{
		for (register int32 j =0; j<A.pixels(); ++j)
		  {
		  fo << A(i,j) << " ";
		  }
		fo << "\n";
		}
	  fo.close();
	  } // END dumpasc
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>




	//***************************************************************
	// * C = dotmult(A,B) = A .* B                                    *
	// * Bert Kampes, 26-Jan-1999                                     *
	// ***************************************************************
	public static <Type> matrix<Type> dotmult(matrix<Type> A, matrix<Type> B)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("dotmult");
	//  #endif
	  matrix<Type> Result = A;
	  Result *= B; // checks are here
	  return Result;
	  } // END dotmult
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = dotdiv(A,B) = A/B                                        *
	// * Bert Kampes, 26-Jan-1999                                     *
	// ***************************************************************
	public static <Type> matrix<Type> dotdiv(matrix<Type> A, matrix<Type> B)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("dotdiv");
	//  #endif
	  matrix<Type> Result = A;
	  Result /= B; // checks are here
	  return Result;
	  } // END dotdiv
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * A = sqr(B)                                                   *
	// * Bert Kampes, 16-Feb-1999                                     *
	// ***************************************************************
	public static <Type> matrix<Type> sqr(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("sqr");
	//  #endif
	  matrix<Type> Result = dotmult(A, A);
	  return Result;
	  } // END sqr
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// *    B = conj(A)                                               *
	// *    Bert Kampes, 02-Mar-1999                                  *
	// ***************************************************************
	public static <Type> matrix<Type> conj(matrix<Type> A)
	  {
	  matrix<Type> Result = A;
	  Result.conj();
	  return Result;
	  } // END conj
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C=diagxmat(vec,B) C=diag(vec) * B;                           *
	// *    Bert Kampes, 22-Feb-1999                                  *
	// ***************************************************************
	public static <Type> matrix<Type> diagxmat(matrix<Type> diag, matrix<Type> B)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("diagxmat");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (min(diag.lines(),diag.pixels()) != 1)
		matERROR.print("diagxmat: sizes A,B: diag is vector.");
	  if (diag.size() != B.lines())
		matERROR.print("diagxmat: sizes A,B: input is vector, matrix.");

	//  #endif
	  matrix<Type> Result =B;
	  if (diag.lines() != 1) // standing
		{
		for (register int32 i =0; i<int32Result.lines(); i++)
		  for (register int32 j =0; j<int32Result.pixels(); j++)
			Result(i,j) *= diag(i,0);
		}
	  else
		{
		for (register int32 i =0; i<int32Result.lines(); i++)
		  for (register int32 j =0; j<int32Result.pixels(); j++)
			Result(i,j) *= diag(0,i);
		}
	  return Result;
	  } // END diagxmat
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * multilook A with factors l,p                                 *
	// * lines(A) have to be multiple of factorL.                     *
	// * size of output is (Al/fl,Ap/fp)                              *
	// *  multilooked is averaged by factorLP.                        *
	// * Bert Kampes, 19-Apr-1999                                     *
	// ***************************************************************
	public static <Type> matrix<Type> multilook(matrix<Type> A, uint factorL, uint factorP)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("multilook.");
		if (A.lines()%factorL) // [MA] we can handle this
		  //matERROR.print("lines A must be multiple of factorL.");
		  matDEBUG.print("For this buffer, lines A is not multiple of factorL.");

	//  #endif
	  if (factorL ==1 && factorP ==1)
		{
		//matrix<Type> R=A;
		//return R;
		return A; // [MA]: fastest solution
		}

	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
		matDEBUG << "multilook input [A] size: " << A.size() << " lines: " << A.lines() << " pixels: " << A.pixels() << " address: " << A; // A.[0] can't be printed if 0
		matDEBUG.print();

	//  #endif
	  if (A.lines()/factorL == 0 || A.pixels()/factorP == 0) // [MA] fix for extra buffer when lines < mlfactor or ...
		{
		DEBUG.print("Multilooking was not necessary for this buffer: buffer.lines() < mlL or buffer.pixels < mlP");
		matrix<Type> R; //=A; // see initialize()
		//R.resize(1,1); // fill with 0 
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//	#if __DEBUGMAT2
		  matDEBUG << "multilook return [R] size: " << R.size() << " lines: " << R.lines() << " pixels: " << R.pixels() << " address: " << R << "\n";
		matDEBUG.print();
	//	#endif
		return R; // NULL
		}

	  Type sum;
	  Type factorLP = new Type(factorL * factorP);
	  //cerr << "multilook: "<< A.lines()/factorL << " " << A.pixels()/factorP << endl;
	  matrix<Type> Result = new matrix(A.lines()/factorL, A.pixels()/factorP);
	  for (register uint i =0; i<Result.lines(); i++)
		{
		for (register uint j =0; j<Result.pixels(); j++)
		  {
		  sum = Type(0.0);
		  for (register uint k =i *factorL; k<(i+1)*factorL; k++)
			{
			for (register uint l =j *factorP; l<(j+1)*factorP; l++)
			  {
			  sum += A(k,l);
			  }
			}
		  Result(i,j) = sum/factorLP;
		  }
		}
	  //cerr << "multilook: l,p " << Result.lines() << " " << Result.pixels() << " size " << Result.size() << endl;
	  return Result;
	  } // END multilook
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * Correlate A with maskB, return C(sizeA)                      *
	// *  egde is set to zero, not pad with zeros                     *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
	public static <Type> matrix<real4> correlate(matrix<Type> A, matrix<Type> Mask)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("Correlate.");
	  matDEBUG.print("not yet correct for complex?: returns real4");
	  if (Mask.lines()<2 || Mask.pixels()<2)
		matERROR.print("very small mask.");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (A.lines()<Mask.lines() || A.pixels()<Mask.pixels())
		matERROR.print("matrix input smaller than mask.");

	//  #endif
	  real8 varM = 0.; // variance of Mask
	  Mask -= mean(Mask);
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntMsk = Mask[0];
	  //for (register uint ii=0; ii<Mask.size(); ii++)
	  //  varM += sqr(*pntMsk++);                         // 1/N later
	  // changed by FvL (for g++/gcc > 4.0):
	  for (register uint ii =0; ii<Mask.size(); ii++)
		{
		varM += sqr(*pntMsk); // 1/N later
		*pntMsk++;
		}

	// ______Compute correlation at these points______
	  uint beginl = (Mask.lines()-1) / 2; // floor
	  uint beginp = (Mask.pixels()-1) / 2; // floor
	  matrix<real4> Result = new matrix(A.lines(), A.pixels()); // init to 0

	// ______First window of A, updated at end of loop______
	  window winA = new window(0, Mask.lines()-1, 0, Mask.pixels()-1);
	  window windef = new window(0,0,0,0); // defaults to total Am

	// ______Correlate part of Result______
	  matrix<Type> Am = new matrix(Mask.lines(), Mask.pixels());
	  for (register uint i =beginl; i<A.lines()-beginl; i++)
		{
		for (register uint j =beginp; j<A.pixels()-beginp; j++)
		  {
		  Am.setdata(windef, A, winA); // Am no allocs.
		  Am -= mean(Am); // center around mean
		  real8 covAM = 0.; // covariance A,Mask
		  real8 varA = 0.; // variance of A(part)
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		  Type *pntM = Mask[0];
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		  Type *pntAm = Am[0];
		  for (register uint l =0; l<Mask.size(); l++)
			{
			//covAM += ((*pntM++) * (*pntAm));        // wait for move pnt
			//varA  += sqr(*pntAm++);                 // pnt ++
			// changed by FvL (for g++/gcc > 4.0):
			covAM += ((*pntM) * (*pntAm)); // wait for move pnt
			varA += sqr(*pntAm); // pnt ++
			*pntM++;
			*pntAm++;
			}
		  // Result(i,j) = covAM / sqrt(varM*varA);
		  Result(i,j) = real4(covAM / Math.sqrt(varM *varA)); // [BO]
		  winA.pixlo++;
		  winA.pixhi++;
		  }
		winA.linelo++;
		winA.linehi++;
		winA.pixlo = 0;
		winA.pixhi = winA.pixlo + Mask.pixels() - 1;
		}
	  return Result;
	  } // END correlate
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * C = -A;                                                      *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
//C++ TO JAVA CONVERTER TODO TASK: Operators cannot be overloaded in Java:
	<Type> matrix<Type> operator -(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("operator -");
	//  #endif
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT1
	  if (!A.size())
		matERROR.print("matrix:: unary minus with empty matrix");
	//  #endif
	  matrix<Type> Result = new matrix(A.lines(), A.pixels());
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntA = A[0];
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntR = Result[0];
	  //for (register uint i=0; i<Result.size(); ++i)
	  //  (*pntR++) = -(*pntA++);
	  // changed by FvL (for g++/gcc > 4.0):
	  for (register uint i =0; i<Result.size(); ++i)
		{
		(*pntR++) = -(*pntA);
		*pntA++;
		}
	  return Result;
	  } // END - (unary)
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * a = mean(A)                                                  *
	// * Bert Kampes, 14-Jan-1999                                     *
	// ***************************************************************
	public static <Type> real8 mean(matrix<Type> A)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("mean");
	//  #endif
	  real8 sum =0.;
	  // ______Ensure stride one memory______
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
	  Type *pntA = A[0];
	  //for (register uint i=0; i<A.size(); ++i)
	  //  sum += (*pntA++);
	  // changed by FvL (for g++/gcc > 4.0):
	  for (register uint i =0; i<A.size(); ++i)
		{
		sum += (*pntA);
		*pntA++;
		}
	  return sum/real8A.size();
	  } // END mean
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * matrix<Type> S = sum(A,dim)                                  *
	// * return sum over dim of A. dim=1 by default.                  *
	// * always returns a matrix, not very handy...                   *
	// * A=[1 2 3]  then: sum(A,1) = [5 7 9]; sum(A,2)=[6]            *
	// *   [4 5 6]                                     [15]           *
	// * Bert Kampes, 28-Mar-2000                                     *
	// ***************************************************************
	public static <Type> matrix<Type> sum(matrix<Type> A, int32 dim)
	  {
	//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
	//  #if __DEBUGMAT2
	  matDEBUG.print("sum");
	//  #endif
	  Type sum = new Type(0);
	  matrix<Type> Res; // may be 1x1 ...
	  if (A.isvector())
		{
		Res.resize(1, 1);
//C++ TO JAVA CONVERTER TODO TASK: Pointer arithmetic is detected on this variable, so pointers on this variable are left unchanged.
		Type *pntA = A[0];
		//for (uint i=0; i<A.size(); i++)
		//  sum += (*pntA++);
		// changed by FvL (for g++/gcc > 4.0):
		for (uint i =0; i<A.size(); i++)
		  {
		  sum += (*pntA);
		  *pntA++;
		  }
		Res(0,0) = sum;
		}
	  else // no vector
		{
		switch (dim)
		  {
		  // ______ sum over rows ______
		  case 1:
			{
			Res.resize(1, A.pixels());
			for (uint i =0; i<A.pixels(); ++i)
			  {
			  for (uint j =0; j<A.lines(); ++j)
				{
				sum += A(j,i);
				}
			  Res(0,i) = sum;
			  sum = Type(0);
			  }
			}
			break;
		  // ______ sum over columns, may be done by pointers for speed ______
		  case 2:
			{
			Res.resize(A.lines(), 1);
			for (uint i =0; i<A.lines(); ++i)
			  {
			  for (uint j =0; j<A.pixels(); ++j)
				{
				sum += A(i,j);
				}
			  Res(i,0) = sum;
			  sum = Type(0);
			  }
			}
			break;
		  default:
			  matERROR.print("sum: dim!={1,2}");
		  }
		} // ifelse vector
	  return Res;
	  } // END sum
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type>



	//***************************************************************
	// * A.mypow(scalar)                                              *
	// * NOmatrix<Type> S = pow(A,scalar)                             *
	// * return pointwize power.                                      *
	// * always returns a matrix, not very handy...                   *
	// #%// BK 26-Oct-2000
	// ***************************************************************
//C++ TO JAVA CONVERTER WARNING: The declaration of the following method implementation was not found:
//ORIGINAL LINE: void matrix<Type>::mypow(Type s)
	public void matrix<Type>mypow(Type s)
	  {
	  for (register uint i =0; i<nrows; ++i)
		for (register uint j =0; j<ncols; ++j)
		  data[i][j] = Math.pow(data[i][j],s);
	  } // END sum
//C++ TO JAVA CONVERTER TODO TASK: The original C++ template specifier was replaced with a Java generic specifier, which may not produce the same behavior:
//ORIGINAL LINE: template <class Type, class Type2>


	//***************************************************************
	// * convert matrix A to type of matrix B                         *
	// * Mahmut Arikan, 07-Jun-2009                                   *
	// * TODO convert_type whole stuff should goto operator = one day.*
	// ***************************************************************
	public static <Type, Type2> void convert_type(matrix<Type> A, matrix<Type2> B)
	  {
	  TRACE_FUNCTION("convert_type(matrix<Type>A-->matrix<Type2>B) (MA 07-Jun-2009)")
	//  #ifdef __DEBUGMAT2
	//    matDEBUG.print("convert_type(A-->B).");
	//  #endif

	  Type pntA = A[0];
	  Type2 pntB = B[0];

	  if (A.lines()!=B.lines() || A.pixels()!=B.pixels())
		{
		DEBUG.print("convert_type aborted since the number of lines or/and pixels of the matrices are not equal.");
		return;
		}
	  else if (A.lines()==0 || A.pixels() == 0)
		{
		DEBUG.print("convert_type aborted since the number of lines or/and pixels of the input matrix is 0.");
		return;
		}

	  //if ( sizeof(Type) == sizeof(Type2) )
	  if (getformat( pntA) == getformat( pntB))
		{
		cerr << "==| convert_type input and output types are the same types" << "==| " << getformat( pntA) << " to " << getformat( pntB) << "\n";
		DEBUG.print("convert_type was not necessary since both input and output types are the same.");
//C++ TO JAVA CONVERTER TODO TASK: The memory management function 'memcpy' has no equivalent in Java:
		memcpy(pntB,pntA,A.size()*sizeof(Type));
		return;
		}

	  for (register uint32 i =0; i<A.size(); i++)
		{
		  pntB = Type2( pntA); // less ambiguous
		   pntB++;
		   pntA++;
		}

	  } // END convert_type
}
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
//***********************************************************************
// * $Source: /users/kampes/DEVELOP/DORIS/doris/src/RCS/matrixbk.cc,v $   *
// * $Revision: 3.14 $                                                    *
// * $Date: 2005/10/18 14:16:32 $                                         *
// * $Author: kampes $                                                    *
// *                                                                      *
// *  template (base) class for matrices.                                 *
// *                                                                      *
// *  ====== Defines (compiler -D option) ======                          *
// *  __USE_VECLIB_LIBRARY__      for using VECLIB library                *
// *  __USE_LAPACK_LIBRARY__      for using LAPACK library                *
// *  __USE_FFTW_LIBRARY__        for using FFTW library                  *
// *  __DEBUGMAT1         index checking, dimensions etc.                 *
// *  __DEBUGMAT2         info, can savely be un-defined                  *
// *                                                                      *
// * Note that data is continuously in memory, because                    *
// *  VECLIB/FFTW assumes this (and me as well...) so a vector of vectors *
// *  cannot be used.                                                     *
// ***********************************************************************


//C++ TO JAVA CONVERTER TODO TASK: There is no preprocessor in Java:
//#if __DEBUGMAT1
//#endif



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
