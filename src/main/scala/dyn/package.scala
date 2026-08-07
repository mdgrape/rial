package rial

/** `dyn` package provides operators whose integer element width is selected at
 *  run time.
 *
 * A module in this package takes a fixed-width input (e.g. 32 bit) and, given a
 * `mode` input, interprets it as a vector of narrower integers (e.g. 8 bit x 4,
 * 16 bit x 2, or a single 32 bit value). The datapath is shared among the modes
 * as much as possible; typically the carry chain is cut at the element boundary
 * that the current mode specifies, so the area is close to that of a single
 * full-width operator.
 */
package object dyn {}
