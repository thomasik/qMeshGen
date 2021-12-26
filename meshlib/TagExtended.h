/////////////////////////////////////////////////////////////////////////////
// TagExtended.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2008-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(TAGEXTENDED_H__INCLUDED)
#define TAGEXTENDED_H__INCLUDED

#include <memory.h>

/**
 * Class defines a linked collection of extra data for mesh points, edges, etc.
 */
class TagExtended
{
public:
	/// Constructor
	TagExtended() : head(0) {}
	/// Copying constructor
	TagExtended(const TagExtended& /* te */) : head(0) { } // -> no copying actually...
	/// Destructor
	~TagExtended(){ removeAllTags(); }
public:
	/// Available types of tags
	enum TagType { TAG_NONE = 0, TAG_ALL, TAG_NEXT,
		TAG_ID, TAG_BID, TAG_MP_2D_3D, TAG_ME_2D_3D, TAG_BOUNDARY_COND, 
		TAG_CUT_MP_2D, TAG_CUT_3D, TAG_CUT_MP_3D, TAG_CUT_MF_3D, TAG_CUT_ME_3D, TAG_MDV, 
		TAG_FIXED, TAG_ACTIVE, TAG_HIDDEN,
		TAG_COLLAPSE_2D, TAG_COLLAPSE_3D, TAG_SPLIT_EDGE_3D,
		TAG_FRONT_1, TAG_FRONT_2,
		TAG_UTC_CUT, TAG_UTC_CUT_PENALTY, TAG_UTC_CUT_DIST, 
		TAG_BOUNDARY_POINT, TAG_BOUNDARY_EDGE, TAG_BOUNDARY_FACE, TAG_SIDE_EDGE, TAG_OUTERHULL_POINT, 
		TAG_MG2D, TAG_TRI_INSPHERE, TAG_MG2D_SM_SWAP, 
		TAG_QUAD, TAG_ACS, TAG_ACS_METRIC,
		TAG_MG3D, TAG_MG3D_SM_SWAP, TAG_MG3D_SM_SWAP_PAR, TAG_MG3D_RECOVERY_DEPTH,
		TAG_DISCRETIZATION_PARAM, TAG_VISUALIZATION, 
		TAG_INSIDE_AREA, TAG_INSIDE_AREA_DIST, TAG_INSIDE_NODE_EDGES, TAG_INSIDE_EDGE_FORBIDDEN, 
		TAG_CUT_ADAPT, TAG_CUT_ORYG_DATA,
		TAG_METRIC_DIFF_RATIO,
		TAG_ADAPT_SURF, TAG_ADAPT_SURF_FIXED, TAG_ADAPT_VOLUME, TAG_ADAPT_CRACK,
		TAG_LAYER_ID, TAG_BLOCK_0, TAG_BLOCK_1, TAG_NORMAL_SP, 
		TAG_LOCAL_SURFACE_DOMAIN, TAG_MAX_CURVATURE, TAG_LOCAL_SURFACE_INNER, TAG_LOCAL_SURFACE_ID,
		TAG_QUALITY, TAG_INV_RANK, TAG_SQQTREE,
//		TAG_LOCAL_SURFACE, TAG_LOCAL_SURFACE_SET, TAG_LOCAL_SURFACE_ORIENT, TAG_LOCAL_CURVE, TAG_LOCAL_CURVE_SET, 
		TAG_GRAIN_NODE_TYPE, TAG_GRAIN_NODE_ID, TAG_GRAIN_NODE_ID2 };
public:
	void copyAllTags(const TagExtended * source) const {
		for(TagData* s = source->head; s; s = s->next){
			bool missing = true;
			for(TagData* t = head; t; t = t->next){
				if(t->data_type == s->data_type){
					t->data = s->data;
					missing = false;
					break;
				}
			}
			if(missing){
				head = new TagData(s->data_type, head);
				head->data = s->data;
			}
		}
	}
	void removeAllTags() const {
		// clear list 
		while(head){
			TagData* temp = head;
			head = head->next;
			delete temp;
		}
		head = 0;
	}
	/// Set/modify (int) tag value
	void setIntTag(TagType tag_type, int value = 1) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type){ t->data.int_value = value; return; }
		}
		// not found -> create
		head = new TagData(tag_type, head);
		head->data.int_value = value;
	}
	/// Set/modify (double) tag value
	void setDoubleTag(TagType tag_type, const double& value) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type){ t->data.double_value = value; return; }
		}
		// not found -> create
		head = new TagData(tag_type, head);
		head->data.double_value = value;
	}
	/// Set/modify (ptr) tag value
	void setPtrTag(TagType tag_type, void* value) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type){ t->data.ptr_value = value; return; }
		}
		// not found -> create
		head = new TagData(tag_type, head);
		head->data.ptr_value = value;
	}
	/// Clear tag -> set zero
	void setZeroTag(TagType tag_type) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type){ t->data.double_value = 0.0; return; }
		}
		// not found -> ok
	}
	/// Whether tag is available
	bool availableTag(TagType tag_type) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return true;
		}
		// not found
		return false;
	}
	/// Whether tag is available and has the given value
	bool checkIntTag(TagType tag_type, int value) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type)
				return t->data.int_value == value;
		}
		// not found
		return false;
	}
	/// Remove tag
	void removeTag(TagType tag_type) const {
		if(!head) return; 
		else if(head->data_type == tag_type){
			TagData* temp = head;
			head = head->next;
			delete temp;
		}else{
			for(TagData* t = head; t->next; t = t->next){
				if(t->next->data_type == tag_type){ 
					TagData* temp = t->next;
					t->next = t->next->next;
					delete temp; 
					return; 
				}
			}
		}
	}
	/// Get (int) tag value
	int getIntTag(TagType tag_type, int default_value = 0) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return t->data.int_value;
		}
		// not found -> return 0
		return default_value;
	}
	/// Get (double) tag value
	double getDoubleTag(TagType tag_type, double default_value = 0.0) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return t->data.double_value;
		}
		// not found -> return 0.0
		return default_value;
	}
	/// Get (ptr) tag value
	void* getPtrTag(TagType tag_type) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return t->data.ptr_value;
		}
		// not found -> return nullptr
		return 0;
	}
	/// Get (tchar) tag value
	bool getTCharTag(TagType tag_type, char tv[8]) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type){
				t->getCharArrayValue(tv);
				return true;
			}
		}
		// not found -> return false
		return false;
	}
	/// Set/modify (tchar) tag value
	void setTCharTag(TagType tag_type,  char tv[8]) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type){ 
				t->setCharArrayValue(tv); 
				return; 
			}
		}
		// not found -> create
		head = new TagData(tag_type, head);
		head->setCharArrayValue(tv);
	}
	/// Set (tchar) flag value
	void setTCharFlag(TagType tag_type,  char tv_i, char tv_f) const {
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type){ 
				t->setCharArrayFlag(tv_i, tv_f); 
				return; 
			}
		}
		// not found -> create
		head = new TagData(tag_type, head);
		head->clearCharArray();
		head->setCharArrayFlag(tv_i, tv_f);
	}
	/// Modifies (int) tag
	void decIntTag(TagType tag_type, int value = 1) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type){ t->data.int_value -= value; return; }
		}
		// not found -> create
		head = new TagData(tag_type, head);
		head->data.int_value = -value;
	}
	/// Modifies (int) tag
	void incIntTag(TagType tag_type, int value = 1) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type){ t->data.int_value += value; return; }
		}
		// not found -> create
		head = new TagData(tag_type, head);
		head->data.int_value = value;
	}
	/// Checks whether the element is marked
	bool nonZeroIntTag(TagType tag_type) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return t->data.int_value != 0;
		}
		return false;
	}
	/// Checks whether the element is not marked
	bool zeroIntTag(TagType tag_type) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return t->data.int_value == 0;
		}
		return true;
	}
	/// Checks whether the element is marked
	bool nonZeroPtrTag(TagType tag_type) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return t->data.ptr_value != 0;
		}
		return false;
	}
	/// Clear flag
	void clearIntFlag(TagType tag_type, int f) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) { t->data.int_value &= ~f; return; }
		}
	}
	/// Set flag
	void setIntFlag(TagType tag_type, int f) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) { t->data.int_value |= f; return; }
		}
		// not found -> create
		head = new TagData(tag_type, head);
		head->data.int_value = f;
	}
	/// Check flag
	bool hasIntFlag(TagType tag_type, int f) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return (t->data.int_value & f) != 0;
		}
		return false;
	}
	/// Check flag
	bool hasAnyIntFlags(TagType tag_type, int f) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return (t->data.int_value & f) != 0;
		}
		return false;
	}
	/// Check flag
	bool hasAllIntFlags(TagType tag_type, int f) const { 
		for(TagData* t = head; t; t = t->next){
			if(t->data_type == tag_type) return (t->data.int_value & f) == f;
		}
		return false;
	}
	///// Printout all tags, for debug
	//void logAllTags() const {
	//	for(TagData* t = head; t; t = t->next){
	//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, t->data_type << " -> " << t->data.int_value << ", ";
	//	}
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;
	//}
protected:
	/// Basic tag node
	class TagData {
	friend class TagExtended;
	public:
		TagData(TagType type, TagData* _next = 0) : data_type(type), next(_next) { }
	public:
		TagType getType() const { return data_type; }
		int getIntValue() const { return data.int_value; }
		double getDoubleValue() const { return data.double_value; }
		void* getPtrValue() const { return data.ptr_value; }
		void getCharArrayValue(char t[8]) const {
			memcpy(t, data.tchar_value, 8);
		}
		void setCharArrayValue(const char t[8]) { 
			memcpy(data.tchar_value, t, 8);
		}
		void clearCharArray() {
			memset(data.tchar_value, 0, 8);
		}
		void setCharArrayFlag(char tv_i, char tv_f) {
			data.tchar_value[tv_i] |= tv_f;
		}
	protected:
		TagType data_type;
		union DataValue {
			int int_value;
			double double_value;
			void* ptr_value;
			char tchar_value[8];
		} data;
		TagData* next;
	};
protected:
	/// simple one-way linked list of tags
	mutable TagData* head;
};

#endif // !defined(TAGEXTENDED_H__INCLUDED)
