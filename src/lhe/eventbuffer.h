#ifndef LHE_EVENTBUFFER_H_
#define LHE_EVENTBUFFER_H_

namespace LHE {

class Event;

class EventBuffer {
  public:
    enum class AllocateStatus {
        SUCCESS,     ///< success
        MALLOCERROR, ///< problem with malloc
        EXISTS       ///< Buffer exists, nothing changed
    };

    enum class AppendStatus {
        SUCCESS,      ///< success
        FULL,         ///< buffer is already full
        MALLOCERROR,  ///< problem with remalloc
        INTERNALERROR ///< the internal stack buffer is too small
    };

    /**
     * @brief constructor
     *
     * Creates an empty buffer. No memory is allocated.
     */
    explicit EventBuffer()
        : data_(NULL), len_(0), cap_(0), n_(0), nmax_(0) {}

    ~EventBuffer() { Free(); }

    /**
     * @brief Allocate Memory
     *
     * Allocate allocates memory for the buffer. NMax is set to 0x7fffffff.
     *
     * @param cap_in_MB    capacity in MB of the buffer
     */
    AllocateStatus Allocate(size_t cap_in_MB);

    /**
     * @brief Append event
     *
     * AppendEvent appends a event to the buffer if the max. number of events is
     * not reached yet.
     * The buffer is resized to fit all events.
     */
    AppendStatus Append(const Event &);

    /**
     * @brief Clear Buffer
     *
     * Clear resets the length of the buffer to 0. Memory is not freed and the
     * actual content is not changed.
     */
    void Clear() {
        len_ = 0;
        n_ = 0;
    }

    /**
     * @brief Free allocated memory
     *
     * Free frees the allocated memory. After a call to Free the object is in
     * the same state as after construction.
     */
    void Free() {
        if(data_ != NULL) {
            free(data_);
            data_ = NULL;
        }
        len_ = 0;
        cap_ = 0;
        n_ = 0;
        nmax_ = 0;
    }

    size_t Length() const { return len_; }
    char *Data() { return data_; }
    const char *Data() const { return data_; }
    size_t N() const { return n_; }
    size_t NMax() const { return nmax_; }

    /**
     * @brief Max number of Events
     *
     * SetNMax set the max. number of events. This function does not remove any
     * events. If there are already more events in this buffer, SetNMax only
     * makes appending impossible.
     */
    void SetNMax(size_t nmax) {
        nmax_ = nmax;
    }
  private:
    char *data_;      ///< data
    size_t len_;      ///< length of the buffer
    size_t grow_cap_; ///< capacity is extended by this if buffer is full
    size_t cap_;      ///< capacity of the buffer
    size_t n_;        ///< number of events in buffer
    size_t nmax_;     ///< max. events
};

} // end namespace LHE

#endif
