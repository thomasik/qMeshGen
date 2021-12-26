/////////////////////////////////////////////////////////////////////////////
// DEquation.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(DEQUATION_H__INCLUDED)
#define DEQUATION_H__INCLUDED

#include "DataHashTable.h"

typedef DataHashTableKeyValue<string, double> DEquationConstTable;

/**
 * This class parses equation from string representation
 *  (creates a linked structure of simple operations)
 *	and evaluetes values and derivatives 
 *	(up to three variables: x, y, z).
 */
class DEquation  
{
public:
	/// Standard constructor
	DEquation(const string& equation_str = "", const DEquationConstTable* ctable = nullptr);
	/// Standard destructor
	virtual ~DEquation();
public:
	/// Available (implemented) operators
	enum DOperator {op_x, op_s = op_x, op_y, op_t = op_y, op_z, op_add, op_sub, op_neg, 
		op_mul, op_div, op_pow, op_sin, op_cos, op_tg, op_PI, op_ln, op_exp, op_sqr, op_sqrt,
		op_not, op_and, op_or, op_lt, op_gt, op_let, op_get, op_eq, op_neq};
	/// Globally used derivatives
	enum DDerivative {deriv_ds, deriv_dt, deriv_du = deriv_ds, deriv_dv = deriv_dt, 
		deriv_dss, deriv_dst, deriv_dtt,
		deriv_dsss, deriv_dsst, deriv_dstt, deriv_dttt};
	/// type of conversion to text
	enum ValueFactor {v_value, v_length, v_area, v_angle, v_angle_rad, v_auto};
public:
	/// Parse string representation of equation, returns true if success
	bool parse(const char* equation_cstr, const DEquationConstTable* ctable = nullptr);
	/// Parse string representation of equation, returns true if success
	bool parse(const string& equation_str, const DEquationConstTable* ctable = nullptr);
	/// Clears the equation structure
	void clear();
	/// Returns the value of the equation for given variables
	double getValue(double x, double y = 0.0, double z = 0.0) const;
	/// Returns a new DEquation object equal to the requested derivative (can be used recursively for derivatives of higher order)
	DEquation* getDerivative(DOperator deriv_op);
	/// Checks, if there are any variables
	bool isConstant() const;
public:
	/// Globally used method of double to string conversion (with suffixes like cm, rad)
	static bool doubleToString(double v, int type, string& text);
	/// Globally used method of double to string conversion (with suffixes like cm, rad)
	static string doubleToString(double v, int type);
	/// Globally used method of string to double conversion (with suffixes like cm, rad), returns true if valid
	static bool stringToDouble(const string& text, int type, double &v);
protected:
	/// Linked list of operations
	struct DOperation{
		/// Operation type
		DOperator op;
		/// Values of operands of this operations (filled by previous operations)
		double data[2];
		/// Reference to operation, where the result of this operation should be placed as operand
		DOperation *result;
		/// Index of operand for result operation
		int	result_id;
		/// Pointers to sub-operations giving the operands for this operation
		DOperation *factors[2];
		/// Working variable
		bool tag;
		/// Next operation to calculate
		DOperation* next;
	} * m_first_operation;
	/// Working stack used during parsing
	struct ParseStack{
		/// Beginning index of substring to parse
		int i1;
		/// Ending index of substring to parse
		int i2;
		/// Reference to result operation
		DOperation* result;
		/// Index of operand in the resulting operation, where the result of this operation should be placed
		int id;
		/// Link to the next substring to process
		ParseStack* next;
	} * m_stack_top;
	/// Working stack used during creation of derivative equations
	struct DerivStack{
		/// Operation to process
		DOperation* operation;
		/// Operation, where result should be placed
		DOperation* result;
		/// Calculation of value or the derivative is requested
		bool deriv;
		/// Index of operand in the resulting operation
		int id;
		/// Link ot the next operation to process
		DerivStack* next;
	} * m_deriv_stack_top;
protected:
	/// Optimizes the equation (removes operation which can be calculated without variables e.g. "0*x", "1+2")
	bool optimize();
	/// Removes spaces from text, returns the count of remaining characters
	size_t clearSpaces(string& text) const;
	/// Returns the type of subexpression in the given substring (e.g. + in "1*x+y")
	int subExpr(string expr, int& from, int& to);
	/// Inserts an operation into the list of operations
	DOperation* insertOperation(DOperator op, DOperation* where, int id);
	/// Takes next substring to parse from stack
	bool pop(int& from, int& to, DOperation *&where, int& id);
	/// Inserts new substring to parse into stack
	void push(int from, int to, DOperation* where, int id);
	/// Takes next operation to "derive" from stack
	bool popDeriv(DOperation* &operation, bool &deriv, DOperation *&where, int& id);
	/// Inserts new operation to "derive" into stack
	void pushDeriv(DOperation* operation, bool deriv, DOperation* where, int id);
protected:
	/// Main operation, where the final result of the whole euqation will be placed
	DOperation m_result;
};

#endif // !defined(DEQUATION_H__INCLUDED)
