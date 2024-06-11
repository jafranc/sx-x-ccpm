/*
 * itf_policies.h
 *
 *  Created on: 17 déc. 2018
 *      Author: jfranc
 */

#ifndef ITF_POLICIES_H_
#define ITF_POLICIES_H_

namespace ccpm{

    /***
     * CPTR template class to implement interfaces
     * @tparam itf_to
     */

//main policy class
    template<class itf_to>
    class interface : public itf_to {

    public:
        interface() : itf_to() {};

    };
}


#endif /* ITF_POLICIES_H_ */
