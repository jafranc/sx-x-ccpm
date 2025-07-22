/*
 * CCPM (Curvature Capture in Porous Media)
 * Copyright (C) 2025 by jacquesn7 (jacquesn7@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 * created by: jacquesn7 (jacquesn7@gmail.com)
 *
 */
#ifndef ITF_POLICIES_H_
#define ITF_POLICIES_H_

namespace ccpm
{

/***
 * CPTR template class to implement interfaces
 * @tparam itf_to
 */

//main policy class
template< class itf_to >
class interface : public itf_to
{

public:
  interface(): itf_to() {};

};
}


#endif /* ITF_POLICIES_H_ */
