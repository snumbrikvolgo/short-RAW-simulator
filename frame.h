#ifndef IITP_FRAME_H
#define IITP_FRAME_H

struct Packet {
    enum Status {
        QUEUED,
        DELIVERED,
        DROPPED_RETRY,
        DROPPED_DELAY
    };

    Packet(double birth_time)
            : birth_time(birth_time), drop_time(-1.0), status(Status::QUEUED)
    { }

    double birth_time;
    double drop_time;
    Status status;
};



#endif //IITP_FRAME_H
