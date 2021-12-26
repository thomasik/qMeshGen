// DEquation.cpp: implementation of the DEquation class.
//////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "common.h"

#include "DEquation.h"

DEquation::DEquation(const string& equation_str, const DEquationConstTable* ctable)
{
	m_first_operation = nullptr;
	m_stack_top = nullptr;
	m_deriv_stack_top = nullptr;
	m_result.data[0] = 0.0;
	if(equation_str != "") parse(equation_str, ctable);
}

DEquation::~DEquation()
{
	clear();
}

void DEquation::clear()
{
	while(m_first_operation){
		DOperation* next = m_first_operation->next;
		delete m_first_operation;
		m_first_operation = next;
	}
	while(m_stack_top){
		ParseStack* next = m_stack_top->next;
		delete m_stack_top;
		m_stack_top = next;
	}
	while(m_deriv_stack_top){
		DerivStack* next = m_deriv_stack_top->next;
		delete m_deriv_stack_top;
		m_deriv_stack_top = next;
	}
}

/// Parse string representation of equation, returns true if success
bool DEquation::parse(const char* equation_cstr, const DEquationConstTable* ctable)
{
	if(equation_cstr == nullptr){
		clear();
		return false;
	}else
		return parse(string(equation_cstr), ctable);
}

bool DEquation::parse(const string& equation_str, const DEquationConstTable* ctable)
{
	const string operators = "+-*/^";
	const string figures = "0123456789.-";

	clear();

	size_t len = equation_str.length();
	if(len < 1 || len > 999) return false;

	string buffer = equation_str;
	for(size_t i = 0; i < len; i++)
		buffer[i] = toupper(buffer[i]);

	len = clearSpaces(buffer);

	int matching = 0;

	for(size_t i = 0; i < len; i++){
		switch(buffer[i]){
		case '(' : ++matching; break;
		case ')' : --matching; break;
		}
		if(matching < 0) return false;
	}
	if(matching > 0) return false;
	if(operators.find(buffer[len-1]) != string::npos) return false;

	DOperation* operation = nullptr;

	push(0, (int)len-1, &m_result, 0);
	int from, to, skip;
	DOperation* result;
	int result_id;
	while(pop(from, to, result, result_id)){
		int i = subExpr(buffer, from, to);
		if(i == -2){	// error
			clear();
			return false;
		}
		if(i == -1 && ctable){ // term -> check additional (optional) table of user-defined constants
			string cname = buffer.substr(from, to-from+1);
			if(ctable->contains(cname)){
				result->data[result_id] = ctable->getValue(cname, 0.0);
				continue;
			}
		}
		switch(i){
		case -1:	// term -> number, function, constant, variable
			switch(buffer[from]){
			case 'X':
			case 'U':
				if(to == from){
					insertOperation(op_x, result, result_id);
				}else{
					clear();
					return false;
				}
				break;
			case 'Y':
			case 'V':
				if(to == from){
					insertOperation(op_y, result, result_id);
				}else{
					clear();
					return false;
				}
				break;
			case 'Z':
				if(to == from){
					insertOperation(op_z, result, result_id);
				}else{
					clear();
					return false;
				}
				break;
			case 'P':
				if(to == from+1 && buffer[to] == 'I'){
					insertOperation(op_PI, result, result_id);
				}else{
					clear();
					return false;
				}
				break;
			case 'S':
				if(buffer.compare(from, 4, "SIN(") == 0){
					operation = insertOperation(op_sin, result, result_id);
					push(from+3, to, operation, 0);
					break;
				}else if(buffer.compare(from, 4, "SQR(") == 0){
					operation = insertOperation(op_sqr, result, result_id);
					push(from+3, to, operation, 0);
					break;
				}else if(buffer.compare(from, 5, "SQRT(") == 0){
					operation = insertOperation(op_sqrt, result, result_id);
					push(from+4, to, operation, 0);
					break;
				}else{
					clear();
					return false;
				}
			case 'C':
				if(buffer.compare(from, 4, "COS(") == 0){
					operation = insertOperation(op_cos, result, result_id);
					push(from+3, to, operation, 0);
					break;
				}else{
					clear();
					return false;
				}
			case 'E':
				if(buffer.compare(from, 4, "EXP(") == 0){
					operation = insertOperation(op_exp, result, result_id);
					push(from+3, to, operation, 0);
					break;
				}else{
					clear();
					return false;
				}
			case 'T':
				if(to == from){
					insertOperation(op_x, result, result_id);	// treat t as a parameter x
					break;
				}else if(buffer.compare(from, 3, "TG(") == 0){
					operation = insertOperation(op_tg, result, result_id);
					push(from+2, to, operation, 0);
					break;
				}else{
					clear();
					return false;
				}
			case 'L':
				if(buffer.compare(from, 3, "LN(") == 0){
					operation = insertOperation(op_ln, result, result_id);
					push(from+2, to, operation, 0);
					break;
				}else{
					clear();
					return false;
				}
			default:
				{
					double d;
					if(stringToDouble(buffer.substr(from, to-from+1), DEquation::v_auto, d)){
						result->data[result_id] = d;
						break;
					}else{
						clear();
						return false;
					}
				}
			}
			break;
		case 0:
			operation = insertOperation(op_neg, result, result_id);
			push(from+1, to, operation, 0);
			break;
		default:
			switch(buffer[i]){
			case '+': operation = insertOperation(op_add, result, result_id); skip = 0; break;
			case '-': operation = insertOperation(op_sub, result, result_id); skip = 0; break;
			case '*': operation = insertOperation(op_mul, result, result_id); skip = 0; break;
			case '/': operation = insertOperation(op_div, result, result_id); skip = 0; break;
			case '^': operation = insertOperation(op_pow, result, result_id); skip = 0; break;
			case '&': operation = insertOperation(op_and, result, result_id); skip = 1; break;
			case '|': operation = insertOperation(op_and, result, result_id); skip = 1; break;
			case '=': operation = insertOperation(op_eq, result, result_id); skip = 1; break;
			case '!': if(buffer[i+1] == '=') {
							operation = insertOperation(op_neq, result, result_id); skip = 1;
					  }else{
							operation = insertOperation(op_not, result, result_id); skip = 0;
					  } break;
			case '>': if(buffer[i+1] == '=') {
							operation = insertOperation(op_get, result, result_id); skip = 1;
					  }else{
							operation = insertOperation(op_gt, result, result_id); skip = 0;
					  } break;
			case '<': if(buffer[i+1] == '=') {
							operation = insertOperation(op_let, result, result_id); skip = 1;
					  }else{
							operation = insertOperation(op_lt, result, result_id); skip = 0;
					  } break;
			default: assert(false);
			}
			if(i > from) push(from, i-1, operation, 0);
			push(i+skip+1, to, operation, 1);
		}
	}

	return optimize();
}

double DEquation::getValue(double x, double y, double z) const
{
	for(DOperation* operation = m_first_operation; operation; operation = operation->next){
		switch(operation->op){
		case op_x:
			operation->result->data[operation->result_id] = x;
			break;
		case op_y:
			operation->result->data[operation->result_id] = y;
			break;
		case op_z:
			operation->result->data[operation->result_id] = z;
			break;
		case op_PI:
			operation->result->data[operation->result_id] = PI;
			break;
		case op_add:
			operation->result->data[operation->result_id] = operation->data[0] + operation->data[1];
			break;
		case op_sub:
			operation->result->data[operation->result_id] = operation->data[0] - operation->data[1];
			break;
		case op_neg:
			operation->result->data[operation->result_id] = -operation->data[0];
			break;
		case op_mul:
			operation->result->data[operation->result_id] = operation->data[0] * operation->data[1];
			break;
		case op_sqr:
			operation->result->data[operation->result_id] = sqr(operation->data[0]);
			break;
		case op_sqrt:
			operation->result->data[operation->result_id] = sqrt(operation->data[0]);
			break;
		case op_div:
			operation->result->data[operation->result_id] = operation->data[0] / operation->data[1];
			break;
		case op_pow:
			if(abs(operation->data[0]) < SMALL_NUMBER){
			operation->result->data[operation->result_id] = 0.0;
			}else{
				operation->result->data[operation->result_id] = exp(operation->data[1] * log(abs(operation->data[0])));
			}
			break;
		case op_sin:
			operation->result->data[operation->result_id] = sin(operation->data[0]);
			break;
		case op_cos:
			operation->result->data[operation->result_id] = cos(operation->data[0]);
			break;
		case op_tg:
			operation->result->data[operation->result_id] = sin(operation->data[0]) / cos(operation->data[0]);
			break;
		case op_exp:
			operation->result->data[operation->result_id] = exp(operation->data[0]);
			break;
		case op_ln:
			operation->result->data[operation->result_id] = log(operation->data[0]);
			break;
		case op_not:
			operation->result->data[operation->result_id] = (operation->data[0] ? 0 : 1);
			break;
		case op_and:
			operation->result->data[operation->result_id] = ((operation->data[0] && operation->data[1]) ? 1 : 0);
			break;
		case op_or:
			operation->result->data[operation->result_id] = ((operation->data[0] || operation->data[1]) ? 1 : 0);
			break;
		case op_gt:
			operation->result->data[operation->result_id] = ((operation->data[0] > operation->data[1]) ? 1 : 0);
			break;
		case op_lt:
			operation->result->data[operation->result_id] = ((operation->data[0] < operation->data[1]) ? 1 : 0);
			break;
		case op_get:
			operation->result->data[operation->result_id] = ((operation->data[0] >= operation->data[1]) ? 1 : 0);
			break;
		case op_let:
			operation->result->data[operation->result_id] = ((operation->data[0] <= operation->data[1]) ? 1 : 0);
			break;
		case op_eq:
			operation->result->data[operation->result_id] = ((operation->data[0] == operation->data[1]) ? 1 : 0);
			break;
		case op_neq:
			operation->result->data[operation->result_id] = ((operation->data[0] == operation->data[1]) ? 0 : 1);
			break;
		default:
			assert(false);
		}
	}

	return m_result.data[0];
}

DEquation* DEquation::getDerivative(DEquation::DOperator deriv_op)
{
	DEquation* derivative = new DEquation();

	if(m_first_operation == nullptr) return derivative;	// equal to 0

	DOperation *operation, *last_operation = nullptr;
	for(operation = m_first_operation; operation; operation = operation->next){
		operation->factors[0] = operation->factors[1] = nullptr;
	}
	for(operation = m_first_operation; operation; operation = operation->next){
		operation->result->factors[operation->result_id] = operation;
		last_operation = operation;
	}

	DOperation *result, *new_operation;
	int result_id;
	bool deriv;

	pushDeriv(last_operation, true, &(derivative->m_result), 0);

	while(popDeriv(operation, deriv, result, result_id)){
		// if deriv !!!!!
		switch(operation->op){
		case op_x:
		case op_y:
		case op_z:
			if(deriv){
				result->data[result_id] = (deriv_op == operation->op) ? 1.0 : 0.0;
			}else{
				derivative->insertOperation(operation->op, result, result_id);
			}
			break;
		case op_PI:
			if(deriv){
				result->data[result_id] = 0.0;
			}else{
				result->data[result_id] = PI;
			}
			break;
		case op_add:	case op_sub:
		case op_not:
		case op_and:	case op_or:
		case op_gt:		case op_lt:
		case op_get:	case op_let:
		case op_eq:		case op_neq:
			if(deriv){
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, new_operation, 0);
				else new_operation->data[0] = 0.0;
				if(operation->factors[1]) pushDeriv(operation->factors[1], true, new_operation, 1);
				else new_operation->data[1] = 0.0;
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
				if(operation->factors[1]) pushDeriv(operation->factors[1], false, new_operation, 1);
				else new_operation->data[1] = operation->data[1];
			}
			break;
		case op_neg:
			if(deriv){
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, new_operation, 0);
				else new_operation->data[0] = 0.0;
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
			}
			break;
		case op_mul:
			if(deriv){
				new_operation = derivative->insertOperation(op_add, result, result_id);
				DOperation* operation_A = derivative->insertOperation(op_mul, new_operation, 0);
				DOperation* operation_B = derivative->insertOperation(op_mul, new_operation, 1);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, operation_A, 0);
				else operation_A->data[0] = 0.0;
				if(operation->factors[1]) pushDeriv(operation->factors[1], false, operation_A, 1);
				else operation_A->data[1] = operation->data[1];
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, operation_B, 0);
				else operation_B->data[0] = operation->data[0];
				if(operation->factors[1]) pushDeriv(operation->factors[1], true, operation_B, 1);
				else operation_B->data[1] = 0.0;
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
				if(operation->factors[1]) pushDeriv(operation->factors[1], false, new_operation, 1);
				else new_operation->data[1] = operation->data[1];
			}
			break;
		case op_sqr:
			if(deriv){
				new_operation = derivative->insertOperation(op_mul, result, result_id);
				new_operation->data[1] = 2.0;
				DOperation* operation_A = derivative->insertOperation(op_mul, new_operation, 0);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, operation_A, 0);
				else operation_A->data[0] = 0.0;
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, operation_A, 1);
				else operation_A->data[1] = operation->data[0];
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
			}
			break;
		case op_sqrt:
			if(deriv){
				if(operation->factors[0]){
					new_operation = derivative->insertOperation(op_div, result, result_id);
					pushDeriv(operation->factors[0], true, new_operation, 0);
					DOperation* operation_A = derivative->insertOperation(op_mul, new_operation, 1);
					operation_A->data[0] = -2.0;
					DOperation* operation_B = derivative->insertOperation(op_sqrt, operation_A, 1);
					pushDeriv(operation->factors[0], false, operation_B, 0);
				}else{
					result->data[result_id] = 0.0;
				}
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
			}
			break;
		case op_div:
			if(deriv){
				new_operation = derivative->insertOperation(op_div, result, result_id);
				DOperation* operation_A = derivative->insertOperation(op_sub, new_operation, 0);
				DOperation* operation_B = derivative->insertOperation(op_sqr, new_operation, 1);
				DOperation* operation_A1 = derivative->insertOperation(op_mul, operation_A, 0);
				DOperation* operation_A2 = derivative->insertOperation(op_mul, operation_A, 1);
				if(operation->factors[1]) pushDeriv(operation->factors[1], false, operation_B, 0);
				else operation_B->data[0] = operation->data[1];
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, operation_A2, 0);
				else operation_A2->data[0] = operation->data[0];
				if(operation->factors[1]) pushDeriv(operation->factors[1], true, operation_A2, 1);
				else operation_A2->data[1] = 0.0;
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, operation_A1, 0);
				else operation_A1->data[0] = 0.0;
				if(operation->factors[1]) pushDeriv(operation->factors[1], false, operation_A1, 1);
				else operation_A1->data[1] = operation->data[1];
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
				if(operation->factors[1]) pushDeriv(operation->factors[1], false, new_operation, 1);
				else new_operation->data[1] = operation->data[1];
			}
			break;
		case op_pow:
			if(deriv){
				new_operation = derivative->insertOperation(op_mul, result, result_id);
				pushDeriv(operation, false, new_operation, 0);
				DOperation* operation_A = derivative->insertOperation(op_add, new_operation, 1);
				DOperation* operation_B = derivative->insertOperation(op_mul, operation_A, 0);
				DOperation* operation_C = derivative->insertOperation(op_mul, operation_A, 1);
				if(operation->factors[1]) pushDeriv(operation->factors[1], false, operation_B, 0);
				else operation_B->data[0] = operation->data[1];
				DOperation* operation_B1 = derivative->insertOperation(op_div, operation_B, 1);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, operation_B1, 0);
				else operation_B1->data[0] = 0.0;
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, operation_B1, 1);
				else operation_B1->data[1] = operation->data[0];
				if(operation->factors[1]) pushDeriv(operation->factors[1], true, operation_C, 0);
				else operation_C->data[0] = 0.0;
				DOperation* operation_C1 = derivative->insertOperation(op_ln, operation_C, 1);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, operation_C1, 0);
				else operation_C1->data[0] = operation->data[0];
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
				if(operation->factors[1]) pushDeriv(operation->factors[1], false, new_operation, 1);
				else new_operation->data[1] = operation->data[1];
			}
			break;
		case op_sin:
			if(deriv){
				new_operation = derivative->insertOperation(op_mul, result, result_id);
				DOperation* operation_A = derivative->insertOperation(op_cos, new_operation, 0);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, new_operation, 1);
				else new_operation->data[1] = 0.0;
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, operation_A, 0);
				else operation_A->data[0] = operation->data[0];
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
			}
			break;
		case op_cos:
			if(deriv){
				new_operation = derivative->insertOperation(op_mul, result, result_id);
				DOperation* operation_A = derivative->insertOperation(op_neg, new_operation, 0);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, new_operation, 1);
				else new_operation->data[1] = 0.0;
				DOperation* operation_A1 = derivative->insertOperation(op_sin, operation_A, 0);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, operation_A1, 0);
				else operation_A1->data[0] = operation->data[0];
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
			}
			break;
		case op_tg:
			if(deriv){
				new_operation = derivative->insertOperation(op_div, result, result_id);
				DOperation* operation_A = derivative->insertOperation(op_neg, new_operation, 0);
				DOperation* operation_B = derivative->insertOperation(op_sqr, new_operation, 1);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, operation_A, 0);
				else operation_A->data[0] = 0.0;
				DOperation* operation_B1 = derivative->insertOperation(op_cos, operation_B, 0);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, operation_B1, 0);
				else operation_B1->data[0] = operation->data[0];
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
			}
			break;
		case op_exp:
			if(deriv){
				new_operation = derivative->insertOperation(op_mul, result, result_id);
				pushDeriv(operation, false, new_operation, 0);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, new_operation, 1);
				else new_operation->data[1] = 0.0;
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
			}
			break;
		case op_ln:
			if(deriv){
				new_operation = derivative->insertOperation(op_div, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], true, new_operation, 0);
				else new_operation->data[0] = 0.0;
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 1);
				else new_operation->data[1] = operation->data[0];
			}else{
				new_operation = derivative->insertOperation(operation->op, result, result_id);
				if(operation->factors[0]) pushDeriv(operation->factors[0], false, new_operation, 0);
				else new_operation->data[0] = operation->data[0];
			}
			break;
		default:
			assert(false);
		}
	}

	if(derivative->optimize()){
		return derivative;
	}else{
		return nullptr;
	}
}

void DEquation::push(int from, int to, DOperation *where, int id)
{
	ParseStack *stack = new ParseStack;
	stack->i1 = from;
	stack->i2 = to;
	stack->result = where;
	stack->id = id;
	stack->next = m_stack_top;
	m_stack_top = stack;
}

bool DEquation::pop(int &from, int &to, DOperation *&where, int &id)
{
	if(m_stack_top == nullptr)
		return false;
	from = m_stack_top->i1;
	to = m_stack_top->i2;
	where = m_stack_top->result;
	id = m_stack_top->id;
	ParseStack* stack = m_stack_top->next;
	delete m_stack_top;
	m_stack_top = stack;
	return true;
}

DEquation::DOperation* DEquation::insertOperation(DOperator op, DOperation *where, int id)
{
	DOperation* operation = new DOperation;
	operation->op = op;
	operation->result = where;
	operation->result_id = id;
	operation->next = m_first_operation;
	return m_first_operation = operation;
}

void DEquation::pushDeriv(DOperation* operation, bool deriv, DOperation* where, int id)
{
	DerivStack *stack = new DerivStack;
	stack->operation = operation;
	stack->deriv = deriv;
	stack->result = where;
	stack->id = id;
	stack->next = m_deriv_stack_top;
	m_deriv_stack_top = stack;
}

bool DEquation::popDeriv(DOperation* &operation, bool &deriv, DOperation *&where, int& id)
{
	if(m_deriv_stack_top == nullptr)
		return false;
	operation = m_deriv_stack_top->operation;
	deriv = m_deriv_stack_top->deriv;
	where = m_deriv_stack_top->result;
	id = m_deriv_stack_top->id;
	DerivStack* stack = m_deriv_stack_top->next;
	delete m_deriv_stack_top;
	m_deriv_stack_top = stack;
	return true;
}

int DEquation::subExpr(string expr, int& from, int& to)
{
	string operators = "+-*/^";
	string figures = "0123456789.-";
	// remove external parenthesis
	int i, j = 0;
	while(expr[from]=='(' && expr[to]==')'){
		for(i = from; i <= to; i++){
			switch(expr[i]){
			case '(' : ++j; break;
			case ')' : --j; break;
			}
			if(j == 0) break;
		}
		if(i == to){
			++from;
			--to;
		}else break;
	}
	int m[7] = { -1, -1, -1, -1, -1, -1, -1 };

	j = 0;
	for(i = from; i <= to; i++){
		switch(expr[i]){
        case '(' : ++j; break;
        case ')' : --j; break;
		}
		if(j == 0){
			switch(expr[i]){
			case '-':
				if(i == from) m[6] = 0;
				else if(i > from+1 && expr[i-1] == 'E' && figures.find(expr[i-2]) != string::npos){
					// nop -> value in exponential form
//					LOG4CPLUS_INFO(MeshLog::logger_console, "exponent", expr.substr(from, to-from+1));
				}else{
					if(operators.find(expr[i-1]) != string::npos) return -2;
					m[3] = i;
				}
				break;
			case '+':
				if(i == from) return -2;
				else{
					if(operators.find(expr[i-1]) != string::npos) return -2;
					m[3] = i;
				}
				break;
			case '*':
			case '/':
				if(i == from) return -2;
				else{
					if(operators.find(expr[i-1]) != string::npos) return -2;
					m[4] = i;
				}
				break;
			case '^':
				if(i == from) return -2;
				else{
					if(operators.find(expr[i-1]) != string::npos) return -2;
					m[5] = i;
				}
				break;
			case '!':
				if(i == to) return -2;
				else if(expr[i+1] == '='){ // equality
					if(i == from || i+1 == to) return -2;
					if(operators.find(expr[i-1]) != string::npos) return -2;
					m[2] = i++;
				}else{ // negation
					if(i == from) m[1] = i;
				}
				break;
			case '&':
			case '|':
				if(i == from || i == to || i+1 == to) return -2;
				else if(expr[i+1] == expr[i]){ // AND -> &&, OR -> ||
					if(operators.find(expr[i-1]) != string::npos) return -2;
					m[0] = i++;
				}else return -2;
				break;
			case '=':
				if(i == from || i == to || i+1 == to) return -2;
				else if(expr[i+1] == expr[i]){ // ==
					if(operators.find(expr[i-1]) != string::npos) return -2;
					m[2] = i++;
				}else return -2;
				break;
			case '>':
			case '<':
				if(i == from || i == to) return -2;
				if(expr[i+1] == '='){ // >=, <=
					if(operators.find(expr[i-1]) != string::npos) return -2;
					if(i+1 == to) return -2;
					else m[2] = i++;
				}else{ // >, <
					if(operators.find(expr[i-1]) != string::npos) return -2;
					m[2] = i;
				}
				break;
			}
		}
	}

	for(i = 0; i < 7; i++)
		if(m[i] >= 0)
			return m[i];

	return -1;
}

string DEquation::doubleToString(double v, int type)
{
	string result;
	DEquation::doubleToString(v, type, result);
	return result;
}

// Zwraca tekstow¹ reprezentacjê wartoœci z uwzglêdnieniem symboli i wspó³czynnika.
// (Wspó³rzêdne podawane w metrach, k¹t w stopniach)
bool DEquation::doubleToString(double v, int type, string& text)
{
	double factors[] = {1.0, 0.01, 1e-3, 1e-6, 1e-9};
	double factors2[] = {1.0, 1e-4, 1e-6, 1e-12, 1e-18};
	const string symbols[] = {"m", "cm", "mm", "mkm", "nm"};
	int index = 0;

	ostringstream str_text;

	switch(type){
	case v_value:
		str_text << v;
		break;
	case v_length:
		if(abs(v) < 0.1){
			while(index < 5 && (abs(v / factors[index]) < 0.1)){
				index++;
			}
			if(index > 4){
				str_text << v;
			}else{
				str_text << (v/factors[index]) << (symbols[index]);
			}
		}else
			str_text << v;
		break;
	case v_area:
		while(index < 5 && (abs(v / factors2[index]) < 0.1)){
			index++;
		}
		if(index > 4) index = 0;
		str_text << (v/factors2[index]) << (symbols[index]) << "2";
		break;
	case v_angle:
		str_text << v << "o";
		break;
	case v_angle_rad:
		str_text << (v * 180.0 / PI) << "o";
		break;
	}

	text = str_text.str();
	return true;
}

// Konwertuje naPIs do liczby: dla wspó³rzêdnych zwraca wartoœæ w metrach, 
//		dla k¹tów - w stopniach.
// Jeœli tekst nie jest prawid³owy zwraca false.
bool DEquation::stringToDouble(const string& text, int type, double &v)
{
	string ending = "";
	istringstream str_text(text);
	
	str_text >> v;
	if(!str_text) return false;
	str_text >> ending;
	//_strlwr(ending);

	switch(type){
	case v_value:	// simple value
		return true;
	case v_length:	// length (position)
		if(ending == ""){
			return true;
		}else if(ending == "m" || ending == "M"){
			return true;
		}else if(ending == "dm" || ending == "DM"){
			v *= 0.1;
		}else if(ending == "cm" || ending == "CM"){
			v *= 0.01;
		}else if(ending == "mm" || ending == "MM"){
			v *= 1e-3;
		}else if(ending == "mkm" || ending == "MKM"){
			v *= 1e-6;
		}else return false;
		break;
	case v_area:	// area
		if(ending == ""){
		}else if(ending == "m2" || ending == "M2"){
		}else if(ending == "dm2" || ending == "DM2"){
			v *= 0.01;
		}else if(ending == "cm2" || ending == "CM2"){
			v *= 1e-4;
		}else if(ending == "mm2" || ending == "MM2"){
			v *= 1e-6;
		}else if(ending == "mkm2" || ending == "MKM2"){
			v *= 1e-12;
		}else return false;
		break;
	case v_angle:	// angle-value
		if(ending == ""){
			// value in degrees -> rad,
			v *= (PI / 180.0);
		}else if(ending == "o"){
			v *= (PI / 180.0);
		}else if(ending == "rad" || ending == "RAD"){
		}else return false;
		break;
	case v_angle_rad:	// angle-value
		break;
	case v_auto:	// ?
		if(ending == ""){
		}else if(ending == "m" || ending == "M"){
		}else if(ending == "dm" || ending == "DM"){
			v *= 0.1;
		}else if(ending == "cm" || ending == "CM"){
			v *= 0.01;
		}else if(ending == "mm" || ending == "MM"){
			v *= 1e-3;
		}else if(ending == "mkm" || ending == "MKM"){
			v *= 1e-6;
		}else if(ending == "m2" || ending == "M2"){
		}else if(ending == "dm2" || ending == "DM2"){
			v *= 0.01;
		}else if(ending == "cm2" || ending == "CM2"){
			v *= 1e-4;
		}else if(ending == "mm2" || ending == "MM2"){
			v *= 1e-6;
		}else if(ending == "mkm2" || ending == "MKM2"){
			v *= 1e-12;
		}else if(ending == "o"){
			v *= (PI / 180.0);
		}else if(ending == "rad" || ending == "RAD"){
		}else return false;
			
		break;
	}
	return true;
}

size_t DEquation::clearSpaces(string& text) const
{
	size_t len = text.length();
	for(size_t i = 0; i < len ; ){
		if(text[i] != ' ')
			++i;
		else{
			text.erase(i,1);
			--len;
		}
	}
	return text.length();
}

bool DEquation::optimize()
{
	DOperation* operation;
	m_result.tag = true;
	for(operation = m_first_operation; operation; operation = operation->next){
		operation->tag = true;
		operation->factors[0] = operation->factors[1] = nullptr;
	}
	for(operation = m_first_operation; operation; operation = operation->next){
		operation->result->factors[operation->result_id] = operation;
	}

	DOperation* last_operation = nullptr;
	operation = m_first_operation;
	while(operation){
		switch(operation->op){
		case op_x: 		case op_y: 		case op_z: 
			operation->tag = operation->result->tag = false;
			break;
		case op_PI:	
			break;
		case op_add:	case op_sub:	case op_mul:
		case op_div:	case op_pow:	case op_sin:
		case op_neg:	case op_cos:	case op_tg:
		case op_exp:	case op_ln:		case op_sqr:
		case op_not:	case op_and:	case op_or:
		case op_gt:		case op_lt:		case op_get:	
		case op_let:	case op_eq:		case op_neq:
		case op_sqrt:
			operation->result->tag &= operation->tag;
			break;
		default:
			assert(false);
		}

		switch(operation->op){
		case op_x:	case op_y:	case op_z:
			break;
		case op_PI:
			operation->result->data[operation->result_id] = PI;
			break;
		case op_add:
			if(operation->tag) 
				operation->result->data[operation->result_id] = operation->data[0] + operation->data[1];
			break;
		case op_sub:
			if(operation->tag)
				operation->result->data[operation->result_id] = operation->data[0]- operation->data[1];
			break;
		case op_neg:
			if(operation->tag)
				operation->result->data[operation->result_id] = - operation->data[0];
			break;
		case op_mul:
			if(operation->tag)
				operation->result->data[operation->result_id] = operation->data[0] * operation->data[1];
			break;
		case op_sqr:
			if(operation->tag)
				operation->result->data[operation->result_id] = sqr(operation->data[0]);
			break;
		case op_sqrt:
			if(operation->tag)
				operation->result->data[operation->result_id] = sqrt(operation->data[0]);
			break;
		case op_not:
			if(operation->tag)
				operation->result->data[operation->result_id] = (operation->data[0] ? 0 : 1);
			break;
		case op_and:
			if(operation->tag)
				operation->result->data[operation->result_id] = ((operation->data[0] && operation->data[1]) ? 1 : 0);
			break;
		case op_or:
			if(operation->tag)
				operation->result->data[operation->result_id] = ((operation->data[0] || operation->data[1]) ? 1 : 0);
			break;
		case op_gt:
			if(operation->tag)
				operation->result->data[operation->result_id] = ((operation->data[0] > operation->data[1]) ? 1 : 0);
			break;
		case op_lt:
			if(operation->tag)
				operation->result->data[operation->result_id] = ((operation->data[0] < operation->data[1]) ? 1 : 0);
			break;
		case op_get:
			if(operation->tag)
				operation->result->data[operation->result_id] = ((operation->data[0] >= operation->data[1]) ? 1 : 0);
			break;
		case op_let:
			if(operation->tag)
				operation->result->data[operation->result_id] = ((operation->data[0] <= operation->data[1]) ? 1 : 0);
			break;
		case op_eq:
			if(operation->tag)
				operation->result->data[operation->result_id] = ((operation->data[0] == operation->data[1]) ? 1 : 0);
			break;
		case op_neq:
			if(operation->tag)
				operation->result->data[operation->result_id] = ((operation->data[0] == operation->data[1]) ? 0 : 1);
			break;
		case op_div:
			if(operation->tag){
				if(operation->data[1] == 0.0){
					clear(); return false;
				}else operation->result->data[operation->result_id] = operation->data[0] / operation->data[1];
			}
			break;
		case op_pow:
			if(operation->tag){				
				if(abs(operation->data[0]) < SMALL_NUMBER){
					operation->result->data[operation->result_id] = 0.0;
				}else{
					operation->result->data[operation->result_id] = exp(operation->data[1] * log(abs(operation->data[0])));
				}
			}
			break;
		case op_sin:
			if(operation->tag)
				operation->result->data[operation->result_id] = sin(operation->data[0]);
			break;
		case op_cos:
			if(operation->tag)
				operation->result->data[operation->result_id] = cos(operation->data[0]);
			break;
		case op_tg:
			if(operation->tag)
				operation->result->data[operation->result_id] = sin(operation->data[0]) / cos(operation->data[0]);
			break;
		case op_exp:
			if(operation->tag)
				operation->result->data[operation->result_id] = exp(operation->data[0]);
			break;
		case op_ln:
			if(operation->tag)
				operation->result->data[operation->result_id] = log(operation->data[0]);
			break;
		default:
			assert(false);
		}

/*	TODO: 0*x, 1*x, 0+x

		if(!operation->tag){
			if(operation->op == op_mul){
				if(
			}
		}
*/
		if(operation->tag){
			if(last_operation) last_operation->next = operation->next;
			else m_first_operation = operation->next;
			delete operation;
			if(last_operation) operation = last_operation->next;
			else if(m_first_operation) operation = m_first_operation;
			else operation = nullptr;
		}else{
			last_operation = operation;
			operation = operation->next;
		}

	}
	return true;
}

bool DEquation::isConstant() const
{
	for(DOperation* operation = m_first_operation; operation; operation = operation->next){
		switch(operation->op){
		case op_x:
		case op_y:
		case op_z:
			return false;
		default:
			break;
		}
	}

	return true;
}
